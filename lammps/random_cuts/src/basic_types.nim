import sequtils

type
  OutputMode* = enum
    Csv = "csv", LmpBalance = "lmp-balance"
  Dimension* = enum
    X, Y, Z
  ZeroOneFloat* = range[0.0..1.0]

converter toZeroOneFloat*(f: Somefloat): ZeroOneFloat =
  f.ZeroOneFloat

converter toZeroOneFloat*(s: seq[Somefloat]): seq[ZeroOneFloat] =
  s.mapIt(it.ZeroOneFloat)


type
  Triple*[T] = tuple
    x, y, z: T

proc prod*[T](t: Triple[T]): T =
  t.x * t.y * t.z

