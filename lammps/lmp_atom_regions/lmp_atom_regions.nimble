# Package

version       = "0.1.0"
author        = "Yannic Kitten"
description   = "lammps-script command generator for custom region cuts and atom creation"
license       = "GPL-3.0-or-later"
srcDir        = "src"
binDir        = "bin"
bin           = @["lmp_atom_regions"]


# Dependencies

requires "nim >= 2.2.4"

task test, "test-examples to execute":
  for cmd in ["bin/lmp_atom_regions tensor 2 2 2 10 10 10 0.5 1000,1,true 0,0.5,1,0,0.5,1,0,0.5,1",
              "bin/lmp_atom_regions staggered 2 2 2 10 10 10 0.0 1000,1,true 0,0.5,1,0,0.25,1,0,0.75,1,0,0.45,1,0,0.33,1,0,0.7,1,0,0.2,1"]:
    echo "cmd: '" & cmd & "'"
    exec cmd
