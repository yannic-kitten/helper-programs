import std / [ strformat, strutils ]

type
  Dimension* = enum
    X, Y, Z
  ZeroOneFloat* = range[0.0..1.0]
  Triple*[T] = tuple
    x, y, z: T

proc parseIntTriple*(args: openArray[string]): Triple[int] {.raises: [ValueError, RangeDefect].} =
  (args[0].parseInt, args[1].parseInt, args[2].parseInt).Triple

proc parseFloatTriple*(args: openArray[string]): Triple[float] {.raises: [ValueError, RangeDefect].} =
  (args[0].parseFloat, args[1].parseFloat, args[2].parseFloat).Triple

proc get*[T](t: Triple[T]; dim: Dimension): T =
  case dim:
    of X: t.x
    of Y: t.y
    of Z: t.z

proc get*[T](t: var Triple[T]; dim: Dimension): var T =
  case dim:
    of X: return t.x
    of Y: return t.y
    of Z: return t.z

proc prod*[T](t: Triple[T]): T =
  t.x * t.y * t.z

proc coord_to_idx*(coord, dim: Triple[int]; mayor: Dimension): int =
  case mayor:
    of X: coord.x * dim.y * dim.z + coord.y * dim.z + coord.z
    of Y: raise newException(ValueError, "Unable to calculate Y-Mayor index, as there is no defined way to do so!")
    of Z: coord.z * dim.y * dim.x + coord.y * dim.x + coord.x



type
  Pair*[T] = tuple
    low, high: T
  RegionBorder* = Triple[Pair[float]]
  Region* = object
    name*: string
    border*: RegionBorder

proc low*(regbrdr: RegionBorder; dim: Dimension): float = regbrdr.get(dim).low
proc low*(regbrdr: var RegionBorder; dim: Dimension): var float = regbrdr.get(dim).low
proc high*(regbrdr: RegionBorder; dim: Dimension): float = regbrdr.get(dim).high
proc high*(regbrdr: var RegionBorder; dim: Dimension): var float = regbrdr.get(dim).high
proc len*(regbrdr: RegionBorder; dim: Dimension): float = regbrdr.high(dim) - regbrdr.low(dim)

proc initRegion*(id: int, name_prefix: string = "reg", border: RegionBorder = ((0.0,1.0), (0.0,1.0), (0.0,1.0))): Region =
  Region(name: name_prefix.strip & fmt"{id:03}", border: border)

proc `*`*(reg: Region; len: Triple[float]): Region =
  result = reg
  result.border.x.low  *= len.x
  result.border.x.high *= len.x
  result.border.y.low  *= len.y
  result.border.y.high *= len.y
  result.border.z.low  *= len.z
  result.border.z.high *= len.z

proc add_gap*(reg: Region; gap: float): Region =
  result = reg
  result.border.x = (reg.border.x.low + gap/2, reg.border.x.high - gap/2)
  result.border.y = (reg.border.y.low + gap/2, reg.border.y.high - gap/2)
  result.border.z = (reg.border.z.low + gap/2, reg.border.z.high - gap/2)



type 
  AtomParams* = tuple
    atoms_per_region: int
    creation_seed: int
    inc_seed: bool

proc parseAtomParams*(args: openArray[string]): AtomParams {.raises: [ValueError].} =
  (args[0].parseInt, args[1].parseInt, args[2].parseBool).AtomParams

