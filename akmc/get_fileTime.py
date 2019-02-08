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
start = int(args[1])
n_jobs = int(args[2])

import os.path, time

state_listdir=[]
for f in glob.glob('*'):
   if os.path.isdir(cwd+'/'+f):
     # state_listdir.append(int(f))
      state_listdir.append(f)
state_listdir.sort()

for f in state_listdir:
 os.chdir(cwd+'/'+str(f))
 for i in range(start, n_jobs):
   #os.system('mkdir '+str(i))
   os.chdir(cwd+'/'+str(f)+'/'+str(i))
   print("last modified: %s" % time.ctime(os.path.getmtime(cwd+'/'+str(f)+'/'+str(i)+'/dynamics.txt')))
