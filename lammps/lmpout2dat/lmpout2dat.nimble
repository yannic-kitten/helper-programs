# Package

version       = "0.1.0"
author        = "Yannic Kitten"
description   = "Converts lammps output to bare data file to be further processed by statistical tools"
license       = "GPL-3.0-or-later"
srcDir        = "src"
binDir        = "bin"
bin           = @["lmpout2dat"]


# Dependencies

requires "nim >= 2.2.4"
