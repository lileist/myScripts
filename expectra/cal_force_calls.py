#!/usr/bin/env python

import sys
import os
import ase,numpy
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
"""
count_lbfgs [start_number] [#ofruns] [rawNumberOfElement]
"""
def main():
    directory = os.getcwd()
    arg = sys.argv
    number = int(arg[2])
    for i in range(int(arg[1]), int(arg[1])+number):
       bh_file=directory+"/run-"+str(i)+'/geo_opt'
       f_bh = open(bh_file,'r')
       output = open(directory+"/run-"+str(i)+'/force_calls', 'w')
       bh_lines = f_bh.readlines()
       force_call = 0
       for line in bh_lines:
          if line.startswith("#"):
             continue
          fields = line.split()
          if fields[1] == '0':
             output.write("%d\n"%(force_call))
             force_call = 1
             continue
          force_call += 1
       print 'run-', i
if __name__ == '__main__':
    main()
    
