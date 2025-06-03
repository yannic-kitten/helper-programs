# Package

version       = "0.2.0"
author        = "Yannic Kitten"
description   = "program to generate random cuts of a system for tensor or staggered grid"
license       = "GPL-3.0-or-later"
srcDir        = "src"
binDir        = "bin"
bin           = @["random_cuts"]


# Dependencies

requires "nim >= 2.2.4"

task release, "Compile with optimization":
  --threads:"off"
  --define:"release"
  setCommand "build"


task examples, "test-examples to execute":
  for cmd in [
              "bin/random_cuts lmp-balance tensor 2 2 2 0.05 78765",
              "bin/random_cuts lmp-balance tensor 4 4 4 0.05 78765",
              "bin/random_cuts csv tensor 2 2 2 0.05 78765",
              "bin/random_cuts csv staggered 2 2 2 0.05 78765",
              "bin/random_cuts csv tiled 2 2 2 0.05 78765",
              ]:
    echo "cmd: '" & cmd & "'"
    exec cmd
    echo ""
