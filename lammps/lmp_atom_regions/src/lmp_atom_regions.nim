import std / [cmdline, strformat, strutils, sequtils, unicode, math]
import general_utils / error
import basic_types, grids

# NOTE: Quickfix arbitrary-cut-order
#   To not spend too much time on integrating an arbitrary cut order into the Grid object structure,
#   all grids remain to be calculated in historic cut order (e.g. zyx for staggered grid).
#   Only after these region borders have benn calculated and can easily accessed by xyz, they are
#   shuffled according to the intended cut order.
#   This is done using the following translation procedure and dimorder attribute of System.
proc changeDimOrder(origBorder: RegionBorder; dimorder: string): RegionBorder =
  if dimorder in ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]:
    for i, d in dimorder:
      let brdr = case i.range[:0..2]:
                  of 0: origBorder.x
                  of 1: origBorder.y
                  of 2: origBorder.z
      case d:
        of 'x': result.x = brdr
        of 'y': result.y = brdr
        of 'z': result.z = brdr
        else: discard
  elif dimorder == "default":
    return origBorder
  else:
    error fmt"invalid dim order '{dimorder}'"

type System = object
  grid: Grid
  gap: float
  gapped_regions: seq[Region]
  atomparams: AtomParams
  dimorder: string  # quickfix (eventually include into grid (as cutorder))

proc initSystem(grid: Grid; gap: float; atomparams: AtomParams; dimorder: string): System =
  result = System(grid: grid, gap: gap, atomparams: atomparams, dimorder: dimorder)
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
    let dimorder = case gridstyle:
      of TensorGridStyle:    args[9].split(',')[0].reversed
      of StaggeredGridStyle: args[9].split(',')[0].reversed
      of TiledGridStyle:     "default"
    let input_cuts = args[9].split(',')[1..^1].map(parseFloat)
    error_if(case gridstyle:
      of TensorGridStyle:    input_cuts.len != (dim.z + 1) + (dim.y + 1) + (dim.x + 1)
      of StaggeredGridStyle: input_cuts.len != (dim.z + 1) + dim.z * (dim.y + 1) + dim.z * dim.y * (dim.x + 1)
      of TiledGridStyle:     input_cuts.len != 2 * (dim.prod() - 1)
      , "The number of cuts does not fit the grid and processor dimensions")
    return initSystem(newGrid(gridstyle, dim, len, input_cuts), region_gap, atomparams, dimorder)
  except:
    error "invalid argument"

proc lmp_region_cmds(sys: System): string =
  for reg in sys.gapped_regions:
    let borders = reg.border.changeDimOrder(sys.dimorder)
    result &= fmt"region {reg.name} block " &
              fmt"{borders.x.low:.6f} {borders.x.high:.6f} " &
              fmt"{borders.y.low:.6f} {borders.y.high:.6f} " &
              fmt"{borders.z.low:.6f} {borders.z.high:.6f}" & "\n"

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
  error.error_info.usage = "./lmp_atom_regions <tensor|staggered|tiled> <nx> <ny> <nz> <lenx> <leny> <lenz> <region-gap> <natoms-per-region,seed,inc-seed?> <dimorder,input-cuts(csv)>"
  main()

