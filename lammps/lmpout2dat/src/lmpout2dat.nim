## Clean up lammps-output-logs in `file.out` and save the actual data (+ header)
##  in `file.dat`, so that it can be read and processed by R as a table of data.
## If the .dat file for a .out file already exists, it is skipped.

import std / [os, osproc, strutils, strformat]

proc dat(file: string): string =
  if file.endswith(".out"):
    return file.replace(".out", ".dat")
  return file & ".dat"

proc main() =
  var files = "ls *.out".execCmdEx.output.strip.splitLines

  for file in files:
    if file.dat.fileExists:
      continue

    var data: seq[string]
    var read = false
    for line in file.lines:
      if (line.strip.startsWith("Step")): read = true   # pre-start of data (header)
      elif line.startsWith("Loop time"): read = false   # post-end of data
      if read:
        data.add line

    if not read:
      file.dat.writeFile(data.join("\n"))
    else: # read should have been set back to false! - Unless an error occoured!
      stderr.writeLine(fmt"Error occoured in: {file}")


when isMainModule:
  main()
