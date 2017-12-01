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
import random

def main():
    arg = sys.argv
    atoms = read_lammps_trj(filename=arg[1], specorder='Au')
    selected_atoms = random.sample(atoms, int(arg[2]))
    wdir = os.getcwd()+'/data/'
    for i in range (0, len(selected_atoms)):
        write_atoms(wdir+str(i)+'.xyz',[selected_atoms[i]])
       
if __name__ == '__main__':
    main()

        
