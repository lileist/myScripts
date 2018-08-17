"""
This code used to get the initial state number, fastest trajectory(start from lowest structure), overall barriers and akmc total time for each transition
How to run:
1. get 'coord_e' first
2. get_transition_stateN.py filenameOFcoord_e max_coord+1
"""
#!/usr/bin/env python

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
   coords_e = dict([int(pair[0]), numpy.array([float(pair[1]), float(pair[2]), float(pair[3])])] for pair in [line.strip().split() for line in f])
with open(state_main_dir+'state_table') as f:
   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])
dynamics = pd.read_table('dynamics.txt', delimiter = r'\s+', skiprows = [0,1], names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])

atoms = None

output=open('e_coord_time.dat','w')
#overall_barrier.write("%10s %6s%3s%-6s %12s %12s %12s %12s\n"%('transition', 'rs','-->', 'ps', 'rs_e', 'ts_e', 'fs_e', 't(akmc)'))
output.write("%14s %12s %12s %4s\n"%('total-time','energy','coords','#ofAu'))
for i in range(len(dynamics.index)):
   try:
      v = coords_e[dynamics['reactant-id'].iloc[i]]
      output.write("%.6E %12.4f %12.4f %4d\n"%(dynamics['total-time'].iloc[i], v[0], v[1], v[2]))
   except:
      break

"""
for i in range(start_coord, end_coord-1):
   akmc_step=dynamics[dynamics['reactant-id']==minima[i-start_coord][0]]['step-number'].iloc[0]
   end=initial[i+1-start_coord][1]

   rs = []
   barrier = []
   filename = current+'/'+str(i)+'TO'+str(i+1)
   output = open(filename+'.dat','w')
   log_structures = Trajectory(filename+'.traj',
                            'w', atoms)
   log_cores = Trajectory(filename+'_cores.traj',
                            'w', atoms)
   selected_dynamics=dynamics[dynamics['step-number']>=akmc_step]
   for j in range(len(selected_dynamics['step-number'])):
      rs_state=selected_dynamics['reactant-id'].iloc[j]
      rs.append(rs_state)
      #print selected_dynamics['step-number'].iloc[j],":", rs
      proc_table=pd.read_table(state_main_dir+str(rs_state)+'/processtable', delimiter = r'\s+', skiprows = [0], names=['proc-id', 'saddle-e', 'prefactor', 'product-id', 'product-e', 'product-prefactor', 'barrier', 'rate', 'repeats'])
      barrier.append(proc_table[proc_table['proc-id']==selected_dynamics['process-id'].iloc[j]]['saddle-e'].iloc[0])
      try:
        index = rs.index(selected_dynamics['product-id'].iloc[j])
        del rs[index:]
        del barrier[index:]
      except:
        pass
      if selected_dynamics['product-id'].iloc[j]==end:
         rs.append(selected_dynamics['product-id'].iloc[j])
         t_akmc = selected_dynamics['total-time'].iloc[j]
         #print "reached", end
         break
   overall_barrier.write("%10s %6d%3s%-6d %12.4f %12.4f %12.4f %.6E\n"%(str(i)+'-->'+str(i+1), rs_state,'-->', end, minima[i-start_coord][1], max(barrier)-minima[i-start_coord][1], states_e[end], Decimal(str(t_akmc))))
   for j in range(len(rs)):
       state_n = rs[j]
       try:
         output.write("%8d  %4.2f %12.4f\n"%(state_n, float(j), states_e[state_n]))
         output.write("%8d  %4.2f %12.4f\n"%(state_n+0.5, float(j+0.5), barrier[j]))
       except:
         pass
       os.chdir(state_main_dir+str(state_n))
       atoms = read('reactant.con',index=0)
       Au_coord, surface_Au, surface_atom = get_coord(atoms)
       log_structures.write(atoms)
       del atoms[surface_atom]
       log_cores.write(atoms)
"""
