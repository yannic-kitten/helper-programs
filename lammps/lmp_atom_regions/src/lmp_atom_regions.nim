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


#[
import std / [sequtils, strutils, strformat, cmdline, sugar, math, deques]
import general_utils / error

type
  GridStyle = enum
    Tensor = "tensor", Staggered = "staggered", Tiled = "tiled"
  IncompatibleGridError = object of Exception
  Triple[T] = tuple
    x, y, z: T
  Dimension = enum
    X, Y, Z
  zero_one = range[0.0..1.0]
  Region = object
    name: string
    border: Triple[tuple[low, high: float]]
  RCBCut = tuple
    cutdim: Dimension
    ratio: zero_one
  RCBCuts = seq[RCBCut]
  CutNode = ref object
    id: int
    cut: RCBCut
    parent: CutNode
    kids: tuple[left: CutNode, right: CutNode]
  CutTree = CutNode
  Cuts = object
    case grid: GridStyle:
      of Tensor, Staggered:
        dimcuts: Triple[seq[float]]
      of Tiled:
        reccuts: RCBCuts
    dim: Triple[int]
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

proc get[T](t: Triple[T]; dim: Dimension): T =
  case dim:
    of X: t.x
    of Y: t.y
    of Z: t.z
proc prod[T](t: Triple[T]): T = t.x * t.y * t.z
proc z_major_idx(coord, dim: Triple[int]): int =
  coord.z * dim.y * dim.x + coord.y * dim.x + coord.x



proc add_gap(reg: Region; gap: float): Region =
  result = reg
  result.border.x = (reg.border.x.low + gap/2, reg.border.x.high - gap/2)
  result.border.y = (reg.border.y.low + gap/2, reg.border.y.high - gap/2)
  result.border.z = (reg.border.z.low + gap/2, reg.border.z.high - gap/2)



proc left(node: CutNode): CutNode = node.kids.left
proc left(node: var CutNode): var CutNode = node.kids.left
proc right(node: CutNode): CutNode = node.kids.right
proc right(node: var CutNode): var CutNode = node.kids.right

proc asCutTree(reccuts: RCBCuts): CutTree =
  result = CutNode(id: 0, cut: reccuts[0], parent: nil, kids: (nil, nil))
  var leafs = [result].toDeque
  for id, cut in reccuts[1..^1]:
    leafs.addLast CutNode(id: id, cut: cut, parent: leafs.peekFirst, kids: (nil, nil))
    if id mod 2 == 1: # left kid
      leafs.peekFirst.left = leafs.peekLast
    else: # right kid
      leafs.peekFirst.right = leafs.peekLast
      discard leafs.popFirst

proc asRCBCuts(tree: CutTree): RCBCuts =
  var nodes = [tree].toDeque
  while nodes.len > 0:
    if not nodes.peekFirst.left.isNil:
      nodes.addLast(nodes.peekFirst.left)
    if not nodes.peekFirst.right.isNil:
      nodes.addLast(nodes.peekFirst.right)
    result.add nodes.popFirst.cut

proc `*`(a, b: RCBCut): RCBCut =
  result = a
  result.ratio *= b.ratio

proc `*`(cut: RCBCut; lens: Triple[float]): RCBCut =
  result = cut
  result.ratio *= lens.get(cut.cutdim)

#  most likely not required - or is it?
#proc `*`(tree: CutTree; lens: Triple[float]): CutTree =
#  result = CutTree(id: tree.id, cut: tree.cut * lens, parent: tree.parent, kids: (nil, nil))
#  if not tree.left.isNil:
#    result.left = tree.left * lens
#  if not tree.right.isNil:
#    result.right = tree.right * lens
#

#  comment in again, when scaling is correct again
#proc find(tree: CutTree; id: int): CutNode =
#  if tree == nil or tree.id == id:
#    return tree
#  let left = tree.left.find(id)
#  if not left.isNil:
#    return left
#  return tree.right.find(id)
#
#proc calc_borders(node: CutNode): Region =
#  if not node.parent.isNil:
#    node.parentcalc_borders
#
#proc getRegion(tree: CutTree; id: int): Region =
#  let node = tree.find(id)
#  discard

  

proc initCuts(grid: GridStyle; dim: Triple[int]; zyx_cuts: seq[float]): Cuts =
  let z_cuts = zyx_cuts[0..<(dim.z+1)]
  var y_cuts, x_cuts: seq[float]
  case grid:
    of Tensor:
      y_cuts = zyx_cuts[z_cuts.len..^1][0..<(dim.y+1)]
      x_cuts = zyx_cuts[(z_cuts.len+y_cuts.len)..^1][0..<(dim.x+1)]
      return Cuts(dim: dim, grid: grid, dimcuts: (x_cuts, y_cuts, z_cuts), zero_one_scale: true)
    of Staggered:
      y_cuts = zyx_cuts[z_cuts.len..^1][0..<(dim.z*(dim.y+1))]
      x_cuts = zyx_cuts[(z_cuts.len+y_cuts.len)..^1][0..<(dim.z*dim.y*(dim.x+1))]
      return Cuts(dim: dim, grid: grid, dimcuts: (x_cuts, y_cuts, z_cuts), zero_one_scale: true)
    of Tiled:
      let reccuts = collect:
        for i in countup(0, zyx_cuts.high, 2):
          (zyx_cuts[i].int.Dimension, zyx_cuts[i+1].zero_one).RcbCut
      return Cuts(dim: dim, grid: grid, reccuts: reccuts, zero_one_scale: true)

proc getRegion(cuts: Cuts; coord: Triple[int]): Region = 
  let dim = cuts.dim
  let idx = z_major_idx(coord, dim)
  result.name = fmt"reg{idx:02}"
  # TODO: adapt for tiled
  case cuts.grid:
    of Tensor:
      result.border.z = (low: cuts.dimcuts.z[coord.z],
                        high: cuts.dimcuts.z[coord.z + 1])
      result.border.y = (low: cuts.dimcuts.y[coord.y],
                        high: cuts.dimcuts.y[coord.y + 1])
      result.border.x = (low: cuts.dimcuts.x[coord.x],
                        high: cuts.dimcuts.x[coord.x + 1])
    of Staggered:
      result.border.z = (low: cuts.dimcuts.z[coord.z],
                        high: cuts.dimcuts.z[coord.z + 1])
      result.border.y = (low: cuts.dimcuts.y[coord.z * (dim.y+1) + coord.y],
                        high: cuts.dimcuts.y[coord.z * (dim.y+1) + coord.y + 1])
      result.border.x = (low: cuts.dimcuts.x[coord.z * dim.y * (dim.x+1) + coord.y * (dim.x+1) + coord.x],
                        high: cuts.dimcuts.x[coord.z * dim.y * (dim.x+1) + coord.y * (dim.x+1) + coord.x + 1])
    of Tiled:
      # using idx as counter for binary tree-nodes (left-to-right first, top-to-bottom)
      raise newException(IncompatibleGridError, "Tiled grid does not support this ")
      result = cuts.reccuts.asCutTree.getRegion(idx)
      discard # TODO

proc getRegion(cuts: Cuts; x, y, z: int): Region = cuts.getRegion((x,y,z))


proc `*`(cuts: Cuts; lens: Triple[float]): Cuts =
  if not cuts.zero_one_scale:
    error "The Cuts are not on a 0-1 scale and can't be rescaled anymore!"
  result = cuts
  case cuts.grid:
    of Tensor, Staggered:
      result.dimcuts.z.applyIt(it * lens.z)
      result.dimcuts.y.applyIt(it * lens.y)
      result.dimcuts.x.applyIt(it * lens.x)
    of Tiled:
      result.reccuts.applyIt(it * lens)
  result.zero_one_scale = false



proc calc_scaled_regions(sys: var System; gap: float) =
  sys.regions = newSeqOfCap[Region](sys.dim.prod())
  let scaled_cuts = sys.cuts * sys.lens
  for z in 0..<sys.dim.z:
    for y in 0..<sys.dim.y:
      for x in 0..<sys.dim.x:
        sys.regions.add scaled_cuts.getRegion(x, y, z).add_gap(gap)

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
      of Tiled: true # TODO
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
  error_info.usage = "./lmp_atom_regions <tensor|staggered|tiled> <nx> <ny> <nz> <lenx> <leny> <lenz> <region-gap> <natoms-per-region,seed,inc-seed> <zyx-cuts(csv)>"
  main()

]#
