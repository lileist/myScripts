#!/usr/bin/env python

import sys
import os
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
"""
Find max and min value for values read from a series of file:
find_max_min start totalNumberOfFile columnNumberOfValue
"""
def main():
    directory = os.getcwd()
    arg = sys.argv
    start = int(arg[1])
    number = int(arg[2])
    column = int(arg[3])
    energy = []
    for i in range(start, start+number):
       f_bh = open(directory+"/run-"+str(i)+'/pot/pot_log.dat', 'r')
       bh_lines = f_bh.readlines()
       for line in bh_lines:
           if line.startswith('#'):
              continue
           fields = line.split()
           energy.append(float(fields[column]))
    max_e=max(energy)
    min_e=min(energy)
    print max_e, min_e
    print max_e-min_e
if __name__ == '__main__':
    main()
    
