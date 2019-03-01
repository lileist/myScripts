#!/usr/bin/env python
"""
This code is used to reorder the pos.con file for AKMC usage in EON
"""

import sys
import os
from ase.io import read,write
from ase.constraints import FixAtoms
import numpy as np

workdir = os.getcwd()
args = sys.argv
rs = read(args[1])
n_Pd = len([atom.index for atom in rs if atom.symbol=='Pd'])
n_Au = len(rs) - n_Pd

inputfile = open(args[1],'r')
output = open(args[1].split('.')[0]+'_re.con','w')
n_line = 0
pd = []
au = []
while True:
   line = inputfile.readline()
   n_line += 1
   if not line:
     break
   if n_line == 8 or n_line == 9:
      output.write('%s %s\n'%(line.split()[1], line.split()[0]))
   if n_line < 8:
      output.write(line)
   if 'Pd' in line:
      pd.append(line)
      inputfile.readline()
      pd.append("Coordinates of Component 1\n")
      for j in range(n_Pd):
        pd.append(inputfile.readline())
   if 'Au' in line:
      au.append(line)
      inputfile.readline()
      au.append("Coordinates of Component 2\n")
      for j in range(n_Au):
        au.append(inputfile.readline())
for line in pd:
   output.write(line)
for line in au:
   output.write(line)
output.close()
os.system('cp '+args[1].split('.')[0]+'_re.con geo_opt/pos.con')
os.chdir(workdir+'/geo_opt')
os.system('~/code/eon/bin/eon')
os.system('mv output '+args[1].split('.')[0])
