## Programs
 - atom_regions: allows calculation of individual regions in tensor or staggered grid; output as lammps-input-script commands for creation of distinct regions and atoms within these
 - random_cuts: generates random cuts (region separators) in a tensor grid to an optional seed, that can be used to manually set the initial balance; useful for observing the convergence speed of balancing methods.
    (going to be enhanced)


## Outdated programs
 - region_calc: calculates the cuts for a tensor grid that equidistantly partitions a system of given dimensions; output as lammps-input-script commands for creation of distinct regions and atoms within these
    (functionally replaced by atom_regions)
