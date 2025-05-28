import std / [ random, sequtils, strutils, math, algorithm ]
import basic_types

#******************** Grid ********************
type
  GridStyle* = enum
    TensorGridStyle = "tensor", StaggeredGridStyle = "staggered", TiledGridStyle = "tiled"
  Grid* = ref object of RootObj
    dim: Triple[int]
    min_dist: float
    rng: Rand

proc newGrid*(grid: GridStyle; dim: Triple[int]; min_dist: float): Grid
method genCuts(grid: var Grid) {.base.} = discard
method csv*(grid: Grid): string {.base.} = discard



#******************** TensorGrid ********************
type TensorDimCuts = seq[ZeroOneFloat]
type TensorGrid* = ref object of Grid
  cuts: Triple[TensorDimCuts]

proc isValid(grid: TensorGrid; dimcuts: TensorDimCuts): bool =
  for i in 0..<dimcuts.high:
    if abs(dimcuts[i] - dimcuts[i.succ]) < grid.min_dist:
      return false
  return dimcuts.len > 0

proc genDimCuts(grid: var TensorGrid; len: int): TensorDimCuts =
  while not grid.isValid(result):
    result = 0.0 & newSeqWith(len.pred, grid.rng.rand(0.0..1.0).round(6)).sorted & 1.0

method genCuts(grid: var TensorGrid) =
  grid.cuts.z = grid.genDimCuts(grid.dim.z)
  grid.cuts.y = grid.genDimCuts(grid.dim.y)
  grid.cuts.x = grid.genDimCuts(grid.dim.x)

method csv*(grid: TensorGrid): string =
  [
    grid.cuts.z.mapIt($it).join(","),
    grid.cuts.y.mapIt($it).join(","),
    grid.cuts.x.mapIt($it).join(","),
  ].join(",")



#******************** StaggeredGrid ********************
type StaggeredDimCuts = seq[ZeroOneFloat]
type StaggeredGrid* = ref object of Grid
  cuts: Triple[seq[StaggeredDimCuts]]

proc isValid(grid: StaggeredGrid; dimcuts: StaggeredDimCuts): bool =
  for i in 0..<dimcuts.high:
    if abs(dimcuts[i] - dimcuts[i.succ]) < grid.min_dist:
      return false
  return dimcuts.len > 0

proc genDimCuts(grid: var StaggeredGrid; len: int): StaggeredDimCuts =
  while not grid.isValid(result):
    result = 0.0 & newSeqWith(len.pred, grid.rng.rand(0.0..1.0).round(6)).sorted & 1.0

method genCuts(grid: var StaggeredGrid) =
  grid.cuts.z = newSeqWith(1, grid.genDimCuts(grid.dim.z))
  grid.cuts.y = newSeqWith(grid.dim.z, grid.genDimCuts(grid.dim.y))
  grid.cuts.x = newSeqWith(grid.dim.z * grid.dim.y, grid.genDimCuts(grid.dim.x))

method csv*(grid: StaggeredGrid): string =
  [
    grid.cuts.z.mapIt($it).join(","),
    grid.cuts.y.mapIt($it).join(","),
    grid.cuts.x.mapIt($it).join(","),
  ].join(",")




#******************** TiledGrid ********************
type TiledCut = tuple
  cutdim: Dimension
  ratio: ZeroOneFloat
type TiledGrid* = ref object of Grid
  cuts: seq[TiledCut]

proc genTiledCut(grid: var TiledGrid): TiledCut =
  ( grid.rng.rand(Dimension),
    grid.rng.rand(grid.min_dist..(1.0-grid.min_dist)).ZeroOneFloat ).TiledCut

method genCuts(grid: var TiledGrid) =
  grid.cuts = newSeqWith(grid.dim.prod.pred, grid.genTiledCut)

method csv*(grid: TiledGrid): string =
  grid.cuts.mapIt($it.cutdim.ord & "," & $it.ratio).join(",")



#************ Grid - implementation ***********
proc newGrid*(grid: GridStyle; dim: Triple[int]; min_dist: float): Grid =
  result = case grid:
    of TensorGridStyle: new TensorGrid
    of StaggeredGridStyle: new StaggeredGrid
    of TiledGridStyle: new TiledGrid
  result.dim = dim
  result.min_dist = min_dist
  result.genCuts


