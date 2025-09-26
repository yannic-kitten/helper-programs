import std / [cmdline, strformat, strutils, sequtils, math]
import general_utils / error
import basic_types, grids

type System = object
  grid: Grid
  gap: float
  gapped_regions: seq[Region]
  atomparams: AtomParams

proc initSystem(grid: Grid; gap: float; atomparams: AtomParams): System =
  result = System(grid: grid, gap: gap, atomparams: atomparams)
  result.gapped_regions = result.grid.getRegions.mapIt(it.add_gap(gap))

proc parseSystem(args: seq[string]): System =
  error_if args.len < 10, "missing arguments"
  error_if args.len > 10, "too many arguments"
  try:
    let gridstyle = parseEnum[GridStyle](args[0])
    let dim = args[1..3].parseIntTriple
    let len = args[4..6].parseFloatTriple
    let region_gap = args[7].parseFloat
    let atomparams = args[8].split(',').parseAtomParams
    let input_cuts = args[9].split(',').map(parseFloat)
    error_if(case gridstyle:
      of TensorGridStyle:    input_cuts.len != (dim.z + 1) + (dim.y + 1) + (dim.x + 1)
      of StaggeredGridStyle: input_cuts.len != (dim.z + 1) + dim.z * (dim.y + 1) + dim.z * dim.y * (dim.x + 1)
      of TiledGridStyle:     input_cuts.len != 2 * (dim.prod() - 1)
      , "The number of cuts does not fit the grid and processor dimensions")
    return initSystem(newGrid(gridstyle, dim, len, input_cuts), region_gap, atomparams)
  except:
    error "invalid argument"

proc lmp_region_cmds(sys: System): string =
  for reg in sys.gapped_regions:
    result &= fmt"region {reg.name} block " &
              fmt"{reg.border.x.low:.6f} {reg.border.x.high:.6f} " &
              fmt"{reg.border.y.low:.6f} {reg.border.y.high:.6f} " &
              fmt"{reg.border.z.low:.6f} {reg.border.z.high:.6f}" & "\n"

proc lmp_atomcreate_cmds(sys: System): string =
  var seed = sys.atomparams.creation_seed
  for reg in sys.gapped_regions:
    result &= fmt"create_atoms 1 random {sys.atomparams.atoms_per_region} {seed} {reg.name}" & "\n"
    if sys.atomparams.inc_seed: seed.inc


proc main() =
  let cliargs = commandLineParams()

  if "-h" in cliargs or "--help" in cliargs or cliargs.len == 0: error.usage_quit()

  let sys = cliargs.parseSystem()
  echo sys.lmp_region_cmds & "\n" & sys.lmp_atomcreate_cmds


when isMainModule:
  error.error_info.usage = "./lmp_atom_regions <tensor|staggered|tiled> <nx> <ny> <nz> <lenx> <leny> <lenz> <region-gap> <natoms-per-region,seed,inc-seed?> <input-cuts(csv)>"
  main()

