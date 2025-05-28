import std / [ sugar, sequtils, deques ]
import basic_types

#******************** Grid ********************
type
  GridStyle* = enum
    TensorGridStyle = "tensor", StaggeredGridStyle = "staggered", TiledGridStyle = "tiled"
  Grid* = ref object of RootObj
    dim: Triple[int]
    len: Triple[float]

proc newGrid*(grid: GridStyle; dim: Triple[int]; len: Triple[float]; input_cuts: seq[float]): Grid
method registerCuts(grid: var Grid; input_cuts: seq[float]) {.base.} = discard
method getRegions*(grid: Grid): seq[Region] {.base.} = discard



#******************** TensorGrid ********************
type TensorDimCuts = seq[ZeroOneFloat]
type TensorGrid* = ref object of Grid
  cuts: Triple[TensorDimCuts]

method registerCuts(grid: var TensorGrid, input_cuts: seq[float] #[ zyx ]#) =
  let input_cuts = input_cuts.mapIt(it.ZeroOneFloat)
  let z_cuts = input_cuts[0..<(grid.dim.z+1)]
  let y_cuts = input_cuts[z_cuts.len..^1][0..<(grid.dim.y+1)]
  let x_cuts = input_cuts[(z_cuts.len+y_cuts.len)..^1][0..<(grid.dim.x+1)]
  grid.cuts = (x_cuts, y_cuts, z_cuts)

proc getRegion(grid: TensorGrid; coord: Triple[int]): Region =
  let dim = grid.dim
  result = initRegion(coord_to_idx(coord, dim, Z))
  result.border.z = (low: grid.cuts.z[coord.z],
                    high: grid.cuts.z[coord.z + 1])
  result.border.y = (low: grid.cuts.y[coord.y],
                    high: grid.cuts.y[coord.y + 1])
  result.border.x = (low: grid.cuts.x[coord.x],
                    high: grid.cuts.x[coord.x + 1])
  return result * grid.len

method getRegions*(grid: TensorGrid): seq[Region] =
  collect:
    for z in 0..<grid.dim.z:
      for y in 0..<grid.dim.y:
        for x in 0..<grid.dim.x:
          grid.getRegion((x,y,z))



#******************** StaggeredGrid ********************
type StaggeredDimCuts = seq[ZeroOneFloat]
type StaggeredGrid* = ref object of Grid
  cuts: Triple[StaggeredDimCuts]

method registerCuts(grid: var StaggeredGrid; input_cuts: seq[float] #[ zyx ]#) =
  let input_cuts = input_cuts.mapIt(it.ZeroOneFloat)
  let z_cuts = input_cuts[0..<(grid.dim.z+1)]
  let y_cuts = input_cuts[z_cuts.len..^1][0..<(grid.dim.z*(grid.dim.y+1))]
  let x_cuts = input_cuts[(z_cuts.len+y_cuts.len)..^1][0..<(grid.dim.z*grid.dim.y*(grid.dim.x+1))]
  grid.cuts = (x_cuts, y_cuts, z_cuts)

proc getRegion(grid: StaggeredGrid; coord: Triple[int]): Region =
  let dim = grid.dim
  result = initRegion(coord_to_idx(coord, dim, Z))
  result.border.z = (low: grid.cuts.z[coord.z],
                    high: grid.cuts.z[coord.z + 1])
  result.border.y = (low: grid.cuts.y[coord.z * (dim.y+1) + coord.y],
                    high: grid.cuts.y[coord.z * (dim.y+1) + coord.y + 1])
  result.border.x = (low: grid.cuts.x[coord.z * dim.y * (dim.x+1) + coord.y * (dim.x+1) + coord.x],
                    high: grid.cuts.x[coord.z * dim.y * (dim.x+1) + coord.y * (dim.x+1) + coord.x + 1])
  return result * grid.len

method getRegions*(grid: StaggeredGrid): seq[Region] =
  collect:
    for z in 0..<grid.dim.z:
      for y in 0..<grid.dim.y:
        for x in 0..<grid.dim.x:
          grid.getRegion((x,y,z))





#************** TiledGrid - RCBCutTree *************
type
  TiledCut = tuple
    cutdim: Dimension
    ratio: ZeroOneFloat
  RCBCutNode = ref object
    cut: TiledCut
    links: tuple[left, right: RCBCutNode]
  RCBCutTree = RCBCutNode

proc left(node: RCBCutNode): RCBCutNode = node.links.left
proc left(node: var RCBCutNode): var RCBCutNode = node.links.left
proc right(node: RCBCutNode): RCBCutNode = node.links.right
proc right(node: var RCBCutNode): var RCBCutNode = node.links.right

converter asRCBCutTree(cuts: seq[TiledCut]): RCBCutTree =
  result = RCBCutNode(cut: cuts[0], links: (nil,nil))
  var leafs = [result].toDeque
  for cut in cuts[1..^1]:
    leafs.addLast RCBCutNode(cut: cut, links: (nil,nil))
    if leafs.peekFirst.left.isNil: # left kid
      leafs.peekFirst.left = leafs.peekLast
    else: # right kid
      leafs.peekFirst.right = leafs.peekLast
      discard leafs.popFirst

proc cutby(regbrdr: RegionBorder; tiledcut: TiledCut): tuple[low, high: RegionBorder] =
  let dim = tiledcut.cutdim
  result = (regbrdr, regbrdr)
  result.low.high(dim) = regbrdr.len(dim) * tiledcut.ratio + regbrdr.low(dim)
  result.high.low(dim) = regbrdr.len(dim) * tiledcut.ratio + regbrdr.low(dim)

proc calc_regionborders(tree: RCBCutTree; border: RegionBorder): seq[RegionBorder] =
  if tree.isNil:
    return @[border]
  let kidsborder = border.cutby(tree.cut)
  result = tree.left.calc_regionborders(kidsborder.low)
  result &= tree.right.calc_regionborders(kidsborder.high)

proc calc_regions(tree: RCBCutTree; gridlen: Triple[float]): seq[Region] =
  let systemborder = (x: (0.0,gridlen.x), y: (0.0,gridlen.y), z: (0.0,gridlen.z))
  let borders = tree.calc_regionborders(systemborder)
  result = collect:
    for i, brdr in borders:
      initRegion(i, border = brdr)

#******************** TiledGrid ********************
type TiledGrid* = ref object of Grid
  cuts: RCBCutTree

method registerCuts(grid: var TiledGrid; input_cuts: seq[float] #[ (dim,ratio), ]#) =
  let cuts = collect:
    for i in countup(0, input_cuts.high, 2):
      (input_cuts[i].int.Dimension, input_cuts[i+1].ZeroOneFloat).TiledCut
  grid.cuts = cuts.asRCBCutTree

method getRegions*(grid: TiledGrid): seq[Region] =
  grid.cuts.calc_regions(grid.len)


#[ # Not required anymore! #
converter asTiledCuts(tree: RCBCutTree): TiledCuts =
  var nodes = [tree].toDeque
  while nodes.len > 0:
    if not nodes.peekFirst.left.isNil:
      nodes.addLast(nodes.peekFirst.left)
    if not nodes.peekFirst.right.isNil:
      nodes.addLast(nodes.peekFirst.right)
    result.add nodes.popFirst.cut
]#





#************ Grid - implementation ***********
proc newGrid*(grid: GridStyle; dim: Triple[int]; len: Triple[float]; input_cuts: seq[float]): Grid =
  result = case grid:
    of TensorGridStyle: new TensorGrid
    of StaggeredGridStyle: new StaggeredGrid
    of TiledGridStyle: new TiledGrid
  result.dim = dim
  result.len = len
  result.registerCuts(input_cuts)

