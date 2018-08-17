#!/usr/bin/env python
"""
This code used to run a batch of MD jobs with lammps
Pre-exiting files required: job.sub data_lammps PdAu.set in.lammps
Example:
run_md_jobs.py [Temperature] [start] [numberOFjobs]
"""

import sys
import os
import re
import glob
import ast
import numpy
from ase.io import read
from ase import Atoms
import random


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



cwd = os.getcwd()
args = sys.argv
temperature = float(args[1])
start = int(args[2])
n_jobs = int(args[3])

f = open('in.md', 'r')
lines = f.readlines()
for i in range(start, n_jobs):
  os.system('mkdir '+str(i))
  f_o = open('./'+str(i)+'/in.md','w')
  rand_n = random.randint(100000,9000000)
  for line in lines:  
    if 'velocity all create' in line:
       f_o.write("velocity all create %8.2f %8d rot yes dist gaussian\n"%(temperature,rand_n))
       continue
#    if 'temp/rescale' in line:
#       f_o.write("fix 2 all temp/rescale 1 %8.2f %8.2f 0.3 1.0\n"%(temperature,temperature))
#       continue
    if 'velocity all scale' in line:
       f_o.write("velocity all scale %8.2f\n"%(temperature))
       continue
    if 'nvt' in line:
       f_o.write("fix 3 all nvt temp %8.2f %8.2f 0.2 tchain 3\n"%(temperature,temperature))
       continue
    f_o.write(line)
  os.system('cp job.sub data_lammps PdAu.set'+' ./'+str(i))
  with cd(str(i)):
       os.system('qsub job.sub')
       print os.getcwd(),'is submitted'
