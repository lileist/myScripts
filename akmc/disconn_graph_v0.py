"""
This code is used to constructure disconnectivity graph using data (/states) from eon_akmc runs
Use 'pele' package from Wales group
accuracy = 1e-2
nlevels  = 40
max_state_n = 1000
emax = 
n_Au =  19
minima_to_draw = 6  19
color = red orange   #blue, orange, green, red, purple, brown, pink, gray, olive, cyan
max_color = blue
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


def compareStru():
#    if min1.cords == min2.cords:
#       return True
    return False

def get_coord(atoms):
    surf_Au = 0
    i = nl('i', atoms,
                     {('Au','Au'):3.3,
                     ('Au','Pd'):3.3,
                     ('Pd','Pd'):3.3
                     })
    coord = numpy.bincount(i)
    index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
    for i in range(len(coord)):
       if i in index_Au:
          if coord[i] < 10:
            surf_Au += 1
    return surf_Au 

paras = readinputs('inputs.ini')
accuracy = float(paras['accuracy'])
nlevels = int(paras['nlevels'])
max_color = paras['max_color']
max_state_n = int(paras['max_state_n'])
min_state_n = int(paras['min_state_n'])
try:
  fig_name = paras['figure_name']
except:
  fig_name = None
try:
  emphasize_state_start = int(paras['emphasize_state_start'])
  emphasize_state_end = int(paras['emphasize_state_end'])
except:
  emphasize_state_start = None
  emphasize_state_end = None
try:
  specified_state = int(paras['specified_state'])
except:
  specified_state = None


draw_symbol = ast.literal_eval(paras['draw_symbol'])
check_structure = ast.literal_eval(paras['check_structure'])
check_e = float(paras['check_e'])
try:
  Emax = float(paras['emax'])
except:
  Emax = None



n_Au = int(paras['n_Au'])

#color used to draw minimum
minima_to_draw = paras['minima_to_draw'].split()
color = paras['color'].split()
symbol = paras['symbol'].split()
Au_seg = {}

for i in range (n_Au+1):
    Au_seg[i] = []

try:
   current = paras['work_dir']
   output_dir = os.getcwd()
except:
   current = os.getcwd()
   output_dir = current

state_listdir=[]

#if(len(states) == 0):
state_main_dir = current+"/states/"
os.chdir(state_main_dir)
db = Database(accuracy=accuracy,compareMinima=compareStru())
#read in states with energy as dictionary
with open('state_table') as f:
   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])

for f in glob.glob('*'):
   try:
      state_listdir.append(int(f))
   except:
      continue
state_listdir.sort()

process_names=['proc', 'saddle energy', 'prefactor', 'productID', 'product energy', 'product prefactor', 'barrier', 'rate', 'repeats']
emphasize = []
specified_min = None
for dir in state_listdir:
   if dir < min_state_n:
      continue
   if dir == max_state_n:
      break
   print "state:", dir
   try:
     rs = db.addMinimum(states_e[dir], [dir])
   except:
     continue
   if dir==min_state_n:
      start_state = rs
   if emphasize_state_start:
      if dir >= emphasize_state_start and dir <= emphasize_state_end:
         emphasize.append(rs)
   if dir == specified_state:
      specified_min = rs

   os.chdir(state_main_dir+str(dir))
   
   atoms = read('reactant.con',index=0)
   Au_seg[get_coord(atoms)].append(rs)
 
   process_table = pd.read_table('processtable', delimiter = r'\s+', skiprows=[0], names=process_names)
   #only forward process stored to avoid duplicated process
   if min_state_n > 0:
       selected_procs = process_table[process_table['productID']>=dir]
   else:
       selected_procs = process_table[process_table['productID']>=0]
   for i in range(len(selected_procs['productID'])):
       ps_ID = selected_procs['productID'].iloc[i]
       #print ps_ID
       if ps_ID >= max_state_n:
          continue
       try:
         ps = db.addMinimum(states_e[ps_ID], [ps_ID], max_n_minima=-1)
         ts = db.addTransitionState(selected_procs['saddle energy'].iloc[i], [dir, ps_ID], rs, ps)
       except:
         continue
print "#of states:",len(db.minima())
#for mini in db.minima():
#   print mini.coords
if Emax is None:
   Emax = -1e20
   for ts in db.transition_states():
      if ts.energy > Emax:
         Emax = ts.energy
   print 'max ts:', Emax
#check the structures with energy larger than check_e
#print paras['check_structure']
if check_structure:
   for ts in db.transition_states():
       if ts.energy > check_e:
          print ts.coords
   for rs in db.minima():
       if ts.energy > check_e:
          print rs.coords
graph = database2graph(db)
dg = DisconnectivityGraph(graph,nlevels=nlevels,Emax=Emax+0.05,node_offset=0)
dg.calculate()
dg.plot()
if min_state_n == 0:
   dg.draw_minima([start_state],c='tab:gray')

if len(emphasize) > 0:
   dg.draw_minima(emphasize,marker='o',c='tab:red')

if specified_min is not None:
   dg.draw_minima([specified_min], marker='8',c='tab:green')
#find max number of atoms on surface
max_key = 0
for key in Au_seg:
   if len(Au_seg[key]) >0 and key > max_key:
      max_key = key
print 'Max # of surface Au:',max_key

if paras['draw_symbol']:
   for i in range(len(minima_to_draw)):
      print minima_to_draw[i],symbol[i],color[i] 
#      try:
      dg.draw_minima(Au_seg[int(minima_to_draw[i])],marker=symbol[i],c='tab:'+color[i])
#      except:
#         continue

   if str(max_key) not in minima_to_draw:
      dg.draw_minima(Au_seg[max_key],marker='<',c='tab:'+max_color)
if fig_name:
   plt.savefig(output_dir+'/'+fig_name)
else:
   plt.show()
