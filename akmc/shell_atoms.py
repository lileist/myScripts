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

def get_coord(atoms):
    surf_Au = 0
    surface_atom = []
    i = nl('i', atoms,
                     {('Au','Au'):3.4,
                     ('Au','Pd'):3.4,
                     ('Pd','Pd'):3.4,
                     ('Pd','Pt'):3.4,
                     ('Au','Pt'):3.4,
                     ('Pt','Pt'):3.4
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

def g_r(atoms, index):
   dist_min = 0.00
   de = 0.02
   gr_numb = 400
   dist_gr = []
   numb=0
   for i in range(400):
       dist_gr.append(0.0)
   for i in index:
       for j in index:
           if i==j:
              continue
           ij_dist = atoms.get_distance(i, j)
           for k in range (gr_numb):
               if ij_dist >= dist_min + float(k)*de and ij_dist < dist_min + float((k+1))*de:
                  dist_gr[k] += 1
                  numb+=1
                  break
   output = open('g_r.dat', 'w')
   bond = 0.0
   prob = 0.0
   for i in range(gr_numb):
       dist_gr[i] = float(dist_gr[i])/numb
       r = dist_min + float(i)*de - de*0.5
       if r > 2.0 and r<3.4:
          bond += r*dist_gr[i]
          prob += dist_gr[i]
       output.write('%15.6f  %15.6f\n' % (r, dist_gr[i]))
   output.write("#average bond length: %15.6f\n"%(bond/prob))
   print "average bond length:",bond/prob

args = sys.argv
atoms = read(args[1],index=0)

Au_coord, surface_Au, surface_atom = get_coord(atoms)
#exclude_atoms = [int(args[2])]
index_Pd = [ atom.index for atom in atoms if atom.symbol=='Pd']
index_Pt = [ atom.index for atom in atoms if atom.symbol=='Pt']
index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
print len(surface_atom)
print "shell_Au: ",len(set(surface_atom).intersection(index_Au))
#print set(surface_atom).intersection(index_Pd)
#g_r(atoms, list(set(surface_atom).intersection(index_Pd)))
g_r(atoms, index_Pt)
#       log_structures.write(atoms)
#       del atoms[surface_atom]
#       log_cores.write(atoms)
