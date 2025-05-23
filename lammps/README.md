## Programs
 - lmp_atom_regions: calculates individual regions in tensor or staggered grid. The output are lammps-input-script commands which define distinct regions and creates a number of atoms for each region.
    Interface: `lmp_atom_regions <tensor|staggered> <nx> <ny> <nz> <lenx> <leny> <lenz> <region-gap> <natoms-per-region,seed,inc-seed> <csv-cuts>`
 - random_cuts: generates random cuts (region separators) in a tensor or staggered grid to an optional seed.It can be used to either manually set the initial balance with a lammps command or pass random cuts to the lmp_atom_regions program. This is useful for observing the convergence speed of balancing methods.
    Interface: `random_cuts <csv|lmp-balance> <tensor|staggered> <nx> <ny> <nz> <min-dist> [<seed>]`
 - lmpout2dat: converts lammps output (.out) files to data (.dat) files that only contain the tabular data including a header with column names. These data files can easily be read and processed by statistical tools.


## Outdated programs
 - region_calc: calculates the cuts for a tensor grid that equidistantly (tensor) partitions a system of given dimensions; output as lammps-input-script commands for creation of distinct regions and atoms within these
    (functionally replaced by lmp_atom_regions)
