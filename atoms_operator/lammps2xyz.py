#!/usr/bin/env python
"""
This code is used to extract geometry and exafs for given energy and s (chi_deviation).
How to use:
   dot_selection  [energy] [s] [pl_steps] [bh_steps] [species_number]
"""
import sys
import os
from expectra.lammps_caller import read_lammps_trj
from expectra.io import write_atoms
def main():
    arg = sys.argv
    atoms = read_lammps_trj(filename=arg[1], specorder='Au')
    write_atoms('au55_trj.xyz',atoms)
if __name__ == '__main__':
    main()

        
