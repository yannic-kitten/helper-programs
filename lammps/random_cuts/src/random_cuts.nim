import std / [random, math, sequtils, algorithm, strutils, strformat, cmdline]

proc error(msg: string) =
  stderr.writeLine(msg)
  QuitFailure.quit

type
  GridStyle = enum
    Tensor = "tensor", Staggered = "staggered"
  OutputKind = enum
    Csv = "csv", LmpBalance = "lmp-balance"
  Triple[T] = tuple
    x, y, z: T
  Cut = seq[float]
  Cuts = object
    grid: GridStyle
    dim: Triple[int]
    min_dist: float
    rng: Rand
    dimcuts: Triple[seq[Cut]]

proc sum[T](t: Triple[T]): T =  t.x + t.y + t.z
proc prod[T](t: Triple[T]): T = t.x * t.y * t.z


proc csv(cut: Cut): string =
  cut.map(`$`).join(",")


proc isValid(cuts: Cuts; cut: Cut): bool =
  for i in 0..<cut.high:
    if abs(cut[i] - cut[i.succ]) < cuts.min_dist:
      return false
  return cut.len > 0

proc generateCut(cuts: var Cuts; n: int): Cut =
  while not cuts.isValid(result):
    result = 0.0 & newSeqWith(n.pred, cuts.rng.rand(0.0..1.0).round(6)).sorted & 1.0

proc calc_cuts(cuts: var Cuts) =
  let nz = 1
  let ny = if cuts.grid == Tensor: 1 else: cuts.dim.z
  let nx = if cuts.grid == Tensor: 1 else: cuts.dim.z * cuts.dim.y
  cuts.dimcuts.z = newSeqWith(nz, cuts.generateCut(cuts.dim.z))
  cuts.dimcuts.y = newSeqWith(ny, cuts.generateCut(cuts.dim.y))
  cuts.dimcuts.x = newSeqWith(nx, cuts.generateCut(cuts.dim.x))

proc initCuts(grid: GridStyle; dim: Triple[int]; min_dist: float; seed: int = 0): Cuts =
  let rng = if seed != 0: initRand(seed) else: initRand()
  result = Cuts(grid: grid, dim: dim, min_dist: min_dist, rng: rng)
  result.calc_cuts

proc parseCuts(args: seq[string]): Cuts =
  var k = 0
  let grid = parseEnum[GridStyle](args[k]) ; k.inc
  let dim = (args[k].parseInt, args[k+1].parseInt, args[k+2].parseInt).Triple ; k.inc 3
  let min_dist = args[k].parseFloat ; k.inc
  var seed = if k < args.len: k.inc ; args[k.pred].parseInt else: 0
  assert k == args.len
  initCuts(grid, dim, min_dist, seed)

proc lmp_balance_cmd(cuts: Cuts): string =
  assert cuts.grid == Tensor
  let cs = (x: cuts.dimcuts.x[0][1..^2].join(" "),
            y: cuts.dimcuts.y[0][1..^2].join(" "),
            z: cuts.dimcuts.z[0][1..^2].join(" "))
  fmt"balance 0.9 x {cs.x} y {cs.y} z {cs.z}"

proc csv(cuts: Cuts): string =
  [cuts.dimcuts.z.map(csv).join(","),
   cuts.dimcuts.y.map(csv).join(","),
   cuts.dimcuts.x.map(csv).join(",")].join(",")


proc main() =
  # ./random_cuts <output-kind> <grid> <procs-x> <procs-y> <procs-z>  <min-dist>  [<seed>]
  let cliargs = commandLineParams()
  let outputKind = parseEnum[OutputKind](cliargs[0])
  let cuts = parseCuts(cliargs[1..^1])

  echo case outputKind:
    of Csv:         cuts.csv
    of LmpBalance:  cuts.lmp_balance_cmd

when isMainModule:
  main()
