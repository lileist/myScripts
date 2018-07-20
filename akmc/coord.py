"""
This code used to calculate coordination number of Au and number of Au atoms on the surface
"""
#!/usr/bin/env python

import sys
import os
import re
import glob
import ast
import pandas as pd
import matplotlib.pyplot as plt
from pele.storage import Database
from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read

def readinputs(filename):
    f=open(filename, 'r')
    parameters = {}
    lines=f.readlines()
    for line in lines:
      if line.startswith('#'):
         continue
      fields = line.split('=')
      parameters[fields[0].strip()]=fields[1].replace("\n","").strip()
    return parameters


def get_coord(atoms):
    surf_Au = 0
    i = nl('i', atoms,
                     {('Au','Au'):3.3,
                     ('Au','Pd'):3.3,
                     ('Pd','Pd'):3.3
                     })
    coord = numpy.bincount(i)
    index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
    Au_coord=0
    for i in range(len(coord)):
       if i in index_Au:
          Au_coord += coord[i]
          if coord[i] < 10:
            surf_Au += 1
    return float(Au_coord)/float(len(index_Au)),surf_Au 


current = os.getcwd()
state_main_dir = current+"/states/"
os.chdir(state_main_dir)
#read in states with energy as dictionary
with open('state_table') as f:
   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])

state_listdir=[]
for f in glob.glob('*'):
   try:
      state_listdir.append(int(f))
   except:
      continue
state_listdir.sort()

for dir in state_listdir:
   os.chdir(state_main_dir+str(dir))
   atoms = read('reactant.con',index=0)
   Au_coord, surface_Au = get_coord(atoms)
   print dir, states_e[dir], Au_coord, surface_Au
