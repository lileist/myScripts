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

current = os.getcwd()
state_main_dir = current+"/states/"
os.chdir(state_main_dir)
#read in states with energy as dictionary

state_listdir=[]
for f in glob.glob('*'):
   try:
      state_listdir.append(int(f))
   except:
      continue
state_listdir.sort()

for dir in state_listdir:
   os.chdir(state_main_dir+str(dir))
   if not os.path.isfile('reactant.con'):
      print dir, 'is empty'
