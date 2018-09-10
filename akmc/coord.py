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
#import matplotlib.pyplot as plt
#from pele.storage import Database
#from pele.landscape import database2graph
#from pele.utils.disconnectivity_graph import DisconnectivityGraph
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read
import collections

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
                     {('Au','Au'):3.2,
                     ('Au','Pd'):3.2,
                     ('Pd','Pd'):3.2
                     })
    coord = numpy.bincount(i)
    index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
    Au_coord=0
    for i in range(len(coord)):
       if i in index_Au:
          Au_coord += coord[i]
          if coord[i] < 11:
            surf_Au += 1
    return float(Au_coord)/float(len(index_Au)),surf_Au 

arg = sys.argv
state_coords = {}
if len(arg) > 1:
   if os.path.isfile(arg[1]):
      coord_file = open(arg[1],'r')
      print "reading coords"
      lines = coord_file.readlines()
      for line in lines:
         if 'state' in line:
            continue
         fields=line.split()
         state_coords[int(fields[0])]=numpy.array([float(fields[1]), float(fields[2]), float(fields[3])])
      coord_file.close()
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
   if dir in state_coords:
      continue
   os.chdir(state_main_dir+str(dir))
   atoms = read('reactant.con',index=0)
   Au_coord, surface_Au = get_coord(atoms)
   state_coords[dir]=numpy.array([states_e[dir], Au_coord, surface_Au])
   print dir, states_e[dir], Au_coord, surface_Au
os.chdir(current)
od_dict = collections.OrderedDict(sorted(state_coords.items()))
output = open('coord_e_all.dat','w')
output.write("%8s %12s %12s %4s\n"%('state','energy','coords','#ofAu'))
for k, v in od_dict.iteritems():
   output.write("%8d %12.4f %12.4f %4d\n"%(k, v[0], v[1], v[2]))
output.close()
