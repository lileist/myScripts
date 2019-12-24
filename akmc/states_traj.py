#!/usr/bin/env python
"""
This code is used to construct trajectory file with states
for help:
states_traj.py --h
"""

import sys
import os
import argparse
import pandas as pd
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.neighborlist import neighbor_list as nl
import numpy

def get_coord(atoms):
    surf_Au = 0
    surface_atom = []
    i = nl('i', atoms,
                     {('Au','Au'):3.4,
                     ('Au','Pd'):3.4,
                     ('Pd','Pd'):3.4
                     })
    coord = numpy.bincount(i)
    index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
    Au_coord=0
    for i in range(len(coord)):
       if coord[i] < 11:
          surface_atom.append(i)
       if i in index_Au:
          Au_coord += coord[i]
          if coord[i] < 10:
            surf_Au += 1
    return float(Au_coord)/float(len(index_Au)),surf_Au,surface_atom

def log_atoms(f, atoms):
    natom = len(atoms)
    f.write("%4d\n"%(natom))
    f.write("\n")
    for atom in atoms:
      f.write("%4s  %12.8f %12.8f %12.8f\n"%(atom.symbol, atom.x, atom.y, atom.z))

current = os.getcwd()
state_main_dir = current+"/states/"

with open(state_main_dir+'state_table') as f:
   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])

atoms = None
parser = argparse.ArgumentParser()
parser.add_argument('--states', type=int, nargs='+', metavar='StateNumbers', 
        default=None,
        help='States that used to constructe trajectories')

parser.add_argument('--start', type=int, metavar='StateNumbers', 
        default=None,
        help='the start point of the trajectory')

parser.add_argument('--end', type=int, metavar='StateNumbers', 
        default=None,
        help='the end point of the trajectory')

parser.add_argument('--akmc_step', type=int, metavar='AKMCstepnumber', 
        default=None,
        help='generate trajectories based on dynamics.txt by assigning a start akmc_step and an end point(product ID)')

args = parser.parse_args()
if args.start is not None and args.end is not None:
   log_structures = Trajectory(str(args.start)+'_'+str(args.end)+'.traj',
                            'w', atoms)
   log_cores = Trajectory(str(args.start)+'_'+str(args.end)+'_cores.traj',
                            'w', atoms)
   for dir in range(args.start, args.end):
       try:
         os.chdir(state_main_dir+str(dir))
       except:
         break
       atoms = read('reactant.con',index=0)
       #Au_coord, surface_Au, surface_atom = get_coord(atoms)
       log_structures.write(atoms)
       #del atoms[surface_atom]
       #log_cores.write(atoms)

if args.states:
   log_structures = Trajectory(str(args.states[0])+'_'+str(args.states[-1])+'.traj',
                            'w', atoms)
   log_cores = Trajectory(str(args.states[0])+'_'+str(args.states[-1])+'_cores.traj',
                            'w', atoms)
   for dir in args.states:
       os.chdir(state_main_dir+str(dir))
       atoms = read('reactant.con',index=0)
       Au_coord, surface_Au, surface_atom = get_coord(atoms)
       log_structures.write(atoms)
       del atoms[surface_atom]
       log_cores.write(atoms)

if args.akmc_step is not None and args.end is not None:
   rs = []
   barrier = []
   ttime = []
   output = open(str(args.akmc_step)+'.dat','w')
   log_structures = Trajectory(str(args.akmc_step)+'.traj',
                            'w', atoms)
   log_cores = Trajectory(str(args.akmc_step)+'_cores.traj',
                            'w', atoms)
   dynamics = pd.read_table('dynamics.txt', delimiter = r'\s+', skiprows = [0,1], names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])
   selected_dynamics=dynamics[dynamics['step-number']>=args.akmc_step]
   #state_coord = pd.read_table('coord_e_all.dat', delimiter = r'\s+', skiprows = [0], names=['state','energy','coords','#ofAu'])
   for i in range(len(selected_dynamics['step-number'])):
      rs.append(selected_dynamics['reactant-id'].iloc[i])
      barrier.append(selected_dynamics['barrier'].iloc[i])
      if i==0:
         ttime.append(0)
      else:
         ttime.append(selected_dynamics['total-time'].iloc[i-1])
      try:
        index = rs.index(selected_dynamics['product-id'].iloc[i])
        del rs[index:]
        del barrier[index:]
        del ttime[index:]
      except:
        pass
      if selected_dynamics['product-id'].iloc[i]==args.end:
         rs.append(selected_dynamics['product-id'].iloc[i])
         ttime.append(selected_dynamics['total-time'].iloc[i])
         break
   log_xyz = open('movie.xyz','w')
   cm_o=None
   for i in range(len(rs)):
       state_n = rs[i]
       try:
         output.write("%8d  %4.2f %12.4f\n"%(state_n, float(i), states_e[state_n]))
         output.write("%8d  %4.2f %12.4f\n"%(state_n+0.5, float(i+0.5), barrier[i]+states_e[state_n]))
       except:
         pass
       os.chdir(state_main_dir+str(state_n))
       atoms = read('reactant.con',index=0)
       Au_coord, surface_Au, surface_atom = get_coord(atoms)
       log_structures.write(atoms)
       if cm_o is None:
         cm_o = atoms.get_center_of_mass()
       cm = atoms.get_center_of_mass()
       atoms.translate(cm_o - cm)
       log_atoms(log_xyz, atoms)
       del atoms[surface_atom]
       log_cores.write(atoms)
       print(state_n, float(i), states_e[state_n], Au_coord, surface_Au, ttime[i])
