import std / [sequtils, strutils, strformat, cmdline]
import general_utils / error

type
  GridStyle = enum
    Tensor = "tensor", Staggered = "staggered"
  Triple[T] = tuple
    x, y, z: T
  Region = object
    name: string
    border: Triple[tuple[low, high: float]]
  Cuts = object
    grid: GridStyle
    dim: Triple[int]
    dimcuts: Triple[seq[float]]
    zero_one_scale: bool
  Lammps_params = tuple
    atoms_per_region: int
    creation_seed: int
    inc_seed: bool
  System = object
    dim: Triple[int]
    lens: Triple[float]
    cuts: Cuts
    regions: seq[Region]
    lmp: Lammps_params

proc prod[T](t: Triple[T]): T = t.x * t.y * t.z
proc z_major_idx(coord, dim: Triple[int]): int =
  coord.z * dim.y * dim.x + coord.y * dim.x + coord.x



proc add_gap(reg: Region; gap: float): Region =
  result = reg
  result.border.x = (reg.border.x.low + gap/2, reg.border.x.high - gap/2)
  result.border.y = (reg.border.y.low + gap/2, reg.border.y.high - gap/2)
  result.border.z = (reg.border.z.low + gap/2, reg.border.z.high - gap/2)



proc initCuts(grid: GridStyle; dim: Triple[int]; zyx_cuts: seq[float]): Cuts =
  let z_cuts = zyx_cuts[0..<(dim.z+1)]
  var y_cuts, x_cuts: seq[float]
  case grid:
    of Tensor:
      y_cuts = zyx_cuts[z_cuts.len..^1][0..<(dim.y+1)]
      x_cuts = zyx_cuts[(z_cuts.len+y_cuts.len)..^1][0..<(dim.x+1)]
    of Staggered:
      y_cuts = zyx_cuts[z_cuts.len..^1][0..<(dim.z*(dim.y+1))]
      x_cuts = zyx_cuts[(z_cuts.len+y_cuts.len)..^1][0..<(dim.z*dim.y*(dim.x+1))]
  Cuts(dim: dim, grid: grid, dimcuts: (x_cuts, y_cuts, z_cuts), zero_one_scale: true)

proc `[]`(cuts: Cuts; coord: Triple[int]): Region = 
  result.name = fmt"reg{z_major_idx(coord, cuts.dim):02}"
  let dim = cuts.dim
  result.border.z = (low: cuts.dimcuts.z[coord.z],
                    high: cuts.dimcuts.z[coord.z + 1])
  case cuts.grid:
    of Tensor:
      result.border.y = (low: cuts.dimcuts.y[coord.y],
                        high: cuts.dimcuts.y[coord.y + 1])
      result.border.x = (low: cuts.dimcuts.x[coord.x],
                        high: cuts.dimcuts.x[coord.x + 1])
    of Staggered:
      result.border.y = (low: cuts.dimcuts.y[coord.z * (dim.y+1) + coord.y],
                        high: cuts.dimcuts.y[coord.z * (dim.y+1) + coord.y + 1])
      result.border.x = (low: cuts.dimcuts.x[coord.z * dim.y * (dim.x+1) + coord.y * (dim.x+1) + coord.x],
                        high: cuts.dimcuts.x[coord.z * dim.y * (dim.x+1) + coord.y * (dim.x+1) + coord.x + 1])

proc `[]`(cuts: Cuts; x, y, z: int): Region = cuts[(x,y,z)]

proc `*`(cuts: Cuts; lens: Triple[float]): Cuts =
  if not cuts.zero_one_scale:
    error "The Cuts are not on a 0-1 scale and can't be rescaled anymore!"
  result = cuts
  result.dimcuts.z.applyIt(it * lens.z)
  result.dimcuts.y.applyIt(it * lens.y)
  result.dimcuts.x.applyIt(it * lens.x)
  result.zero_one_scale = false



proc calc_scaled_regions(sys: var System; gap: float) =
  sys.regions = newSeqOfCap[Region](sys.dim.prod())
  let scaled_cuts = sys.cuts * sys.lens
  for z in 0..<sys.dim.z:
    for y in 0..<sys.dim.y:
      for x in 0..<sys.dim.x:
        sys.regions.add scaled_cuts[x, y, z].add_gap(gap)

proc initSystem(dim: Triple[int]; lens: Triple[float]; cuts: Cuts; region_gap = 0.0; lmp: Lammps_params = (1000, 12345, true)): System =
  result = System(dim: dim, lens: lens, cuts: cuts, lmp: lmp)
  result.calc_scaled_regions(region_gap)

proc parseSystem(args: seq[string]): System =
  error_if args.len < 10, "missing arguments"
  error_if args.len > 10, "too many arguments"
  try:
    let grid = parseEnum[GridStyle](args[0])
    let dim = (args[1].parseInt, args[2].parseInt, args[3].parseInt).Triple
    let lens = (args[4].parseFloat, args[5].parseFloat, args[6].parseFloat)
    let region_gap = args[7].parseFloat
    let lmp_raw = args[8].split(',')
    let lmp = (lmp_raw[0].parseInt, lmp_raw[1].parseInt, lmp_raw[2].parseBool)
    let zyx_cuts = args[9].split(',').map(parseFloat)
    error_if(case grid:
      of Tensor:    zyx_cuts.len != (dim.z + 1) + (dim.y + 1) + (dim.x + 1)
      of Staggered: zyx_cuts.len != (dim.z + 1) + dim.z * (dim.y + 1) + dim.z * dim.y * (dim.x + 1)
      , "The number of cuts does not fit the grid and processor dimensions")
    let cuts = initCuts(grid, dim, zyx_cuts)
    return initSystem(dim, lens, cuts, region_gap, lmp)
  except:
    error "invalid argument"

proc lmp_region_cmds(sys: System): string =
  for reg in sys.regions:
    result &= fmt"region {reg.name} block " &
              fmt"{reg.border.x.low:.6f} {reg.border.x.high:.6f} " &
              fmt"{reg.border.y.low:.6f} {reg.border.y.high:.6f} " &
              fmt"{reg.border.z.low:.6f} {reg.border.z.high:.6f}" & "\n"

proc lmp_atomcreate_cmds(sys: System): string =
  var seed = sys.lmp.creation_seed
  for reg in sys.regions:
    result &= fmt"create_atoms 1 random {sys.lmp.atoms_per_region} {seed} {reg.name}" & "\n"
    if sys.lmp.inc_seed: seed.inc



proc main() =
  let cliargs = commandLineParams()

  if "-h" in cliargs or "--help" in cliargs or cliargs.len == 0: usage_quit()

  let sys = parseSystem(cliargs)
  echo sys.lmp_region_cmds & "\n" & sys.lmp_atomcreate_cmds


when isMainModule:
  error_info.usage = "./lmp_atom_regions <tensor|staggered> <nx> <ny> <nz> <lenx> <leny> <lenz> <region-gap> <natoms-per-region,seed,inc-seed> <zyx-cuts(csv)>"
  main()

