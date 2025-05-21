import std / [random, math, sequtils, algorithm, strutils, cmdline]
import general_utils / error

type
  GridStyle = enum
    Tensor = "tensor", Staggered = "staggered"
  OutputMode = enum
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
  let ny = if cuts.grid == Tensor: 1 else: #[Staggered]# cuts.dim.z
  let nx = if cuts.grid == Tensor: 1 else: #[Staggered]# cuts.dim.z * cuts.dim.y
  cuts.dimcuts.z = newSeqWith(nz, cuts.generateCut(cuts.dim.z))
  cuts.dimcuts.y = newSeqWith(ny, cuts.generateCut(cuts.dim.y))
  cuts.dimcuts.x = newSeqWith(nx, cuts.generateCut(cuts.dim.x))

proc initCuts(grid: GridStyle; dim: Triple[int]; min_dist: float; seed: int = 0): Cuts =
  let rng = if seed != 0: initRand(seed) else: initRand()
  result = Cuts(grid: grid, dim: dim, min_dist: min_dist, rng: rng)
  result.calc_cuts

proc parseCuts(args: seq[string]): Cuts =
  error_if args.len < 5, "missing arguments"
  error_if args.len > 6, "too many arguments"
  try:
    let grid = parseEnum[GridStyle](args[0])
    let dim = (args[1].parseInt, args[2].parseInt, args[3].parseInt).Triple
    error_if([dim.x, dim.y, dim.z].anyIt(it == 0), "Only 3D grids are supported! (minimal dims: 1x1x1)")
    let min_dist = args[4].parseFloat
    error_if([1/dim.x, 1/dim.y, 1/dim.z].anyIt(it <= min_dist), "Minimal distance is too large")
    var seed = if args.len > 5: args[5].parseInt else: 0
    return initCuts(grid, dim, min_dist, seed)
  except:
    error "invalid argument"

proc lmp_balance_cmd(cuts: Cuts): string =
  error_if cuts.grid != Tensor, "The `lmp-balance` output mode only works for tensor grid"
  let cs = (x: cuts.dimcuts.x[0][1..^2].join(" "),
            y: cuts.dimcuts.y[0][1..^2].join(" "),
            z: cuts.dimcuts.z[0][1..^2].join(" "))
  if cs.x.len + cs.y.len + cs.z.len == 0:
    return ""
  "balance 0.9 " & (if cs.x.len > 0: " x " & cs.x else: "") &
                   (if cs.y.len > 0: " y " & cs.y else: "") &
                   (if cs.z.len > 0: " z " & cs.z else: "")

proc csv(cuts: Cuts): string =
  [cuts.dimcuts.z.map(csv).join(","),
   cuts.dimcuts.y.map(csv).join(","),
   cuts.dimcuts.x.map(csv).join(",")].join(",")


proc main() =
  let cliargs = commandLineParams()

  if "-h" in cliargs or "--help" in cliargs or cliargs.len == 0: usage_quit()

  var outputMode: OutputMode
  try: outputMode = parseEnum[OutputMode](cliargs[0])
  except: error "invalid argument"

  let cuts = parseCuts(cliargs[1..^1])

  echo case outputMode:
    of Csv:         cuts.csv
    of LmpBalance:  cuts.lmp_balance_cmd


when isMainModule:
  error_info.usage = "./random_cuts <csv|lmp-balance> <tensor|staggered> <nx> <ny> <nz> <min-dist> [<seed>]"
  main()

