#!/usr/bin/env python
"""
This code is used to get the evolution of number of surface Au atoms with time basined on 'dynamics.txt'. 
How to run:
1. get 'coord_e' first
2. coord_evolve.py fileNameOfcoordinates
"""

import sys
import os
import re
import glob
import ast
import pandas as pd
from decimal import Decimal
from pele.storage import Database
from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read, write
from ase.io.trajectory import Trajectory


args = sys.argv
current = os.getcwd()
state_main_dir = current+'/states/'

#read in coordination info
#coord_e = pd.read_table(args[1], delimiter = r'\s+', names=['state-number', 'energy', 'coord-number', 'surface-Au'])
coords_e={}
with open(args[1]) as f:
#   for line in f:
#      fields = line.strip().split()
#      coords_e[int()]
   coords_e = dict([int(pair[0]), numpy.array([float(pair[1]), float(pair[2]), float(pair[3])])] for pair in [line.strip().split() for line in f if 'state' not in line])
with open(state_main_dir+'state_table') as f:
   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])
dynamics = pd.read_table('dynamics.txt', delimiter = r'\s+', skiprows = [0,1], names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])

atoms = None

output=open('e_coord_time.dat','w')
#overall_barrier.write("%10s %6s%3s%-6s %12s %12s %12s %12s\n"%('transition', 'rs','-->', 'ps', 'rs_e', 'ts_e', 'fs_e', 't(akmc)'))
output.write("%14s %12s %12s %4s\n"%('total-time','energy','coords','#ofAu'))
for i in range(len(dynamics.index)):
   try:
      v_r = coords_e[dynamics['reactant-id'].iloc[i]]
      v_p = coords_e[dynamics['product-id'].iloc[i]]
      output.write("%.6E %12.4f %12.4f %4d\n"%(dynamics['total-time'].iloc[i], v_r[0], v_r[1], v_r[2]))
      output.write("%.6E %12.4f %12.4f %4d\n"%(dynamics['total-time'].iloc[i], v_p[0], v_p[1], v_p[2]))
   except:
      break

