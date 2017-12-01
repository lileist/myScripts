#!/usr/bin/env python
from ase.io import read, write
from expectra.lammps_caller import read_lammps_trj
import sys


def main():
   arg = sys.argv
   inputfile = arg[1]
   spec_one = arg[2]
   print [spec_one]
   spec_two = arg[3]
   p1 = read_lammps_trj(filename=inputfile,
                        skip = 1,
                        every = 1,
                        specorder = [spec_one, spec_two])
   write(filename='au.xyz',images=p1,format='xyz')
if __name__ == '__main__':
    main()
