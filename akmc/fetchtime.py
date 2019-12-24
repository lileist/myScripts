#!/usr/bin/env python
"""
This code is used to get the evolution of number of surface Au atoms with time basined on 'dynamics.txt'. 
How to run:
1. get 'coord_e' first
2. coord_evolve.py fileNameOfcoordinates
"""

import sys
import os
import pandas as pd
import numpy



#read in coordination info
#coord_e = pd.read_table(args[1], delimiter = r'\s+', names=['state-number', 'energy', 'coord-number', 'surface-Au'])
dynamics = pd.read_table('dynamics.txt', delimiter = r'\s+', skiprows = [0,1], names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])

atoms = None

output=open('e_coord_time.dat','w')
#overall_barrier.write("%10s %6s%3s%-6s %12s %12s %12s %12s\n"%('transition', 'rs','-->', 'ps', 'rs_e', 'ts_e', 'fs_e', 't(akmc)'))
output.write("%14s %12s %12s %4s\n"%('total-time','energy','coords','#ofAu'))
states = numpy.loadtxt('./states/state_dft')[:,1]
for i in range(len(dynamics.index)):
   try:
      #output.write("%.6E %12.4f \n"%(dynamics['total-time'].iloc[i]*10**12, dynamics['energy'].iloc[i]))
      #output.write("%.6E %12.4f \n"%(dynamics['total-time'].iloc[i+1]*10**12,dynamics['energy'].iloc[i]))
      output.write("%.6E %12.4f %12.4f\n"%(dynamics['total-time'].iloc[i]*10**12, states[int(dynamics['reactant-id'].iloc[i])], dynamics['energy'].iloc[i]))
      output.write("%.6E %12.4f %12.4f\n"%(dynamics['total-time'].iloc[i+1]*10**12,states[int(dynamics['reactant-id'].iloc[i])], dynamics['energy'].iloc[i]))
   except:
      break

