"""
This code used to get the initial state number, fastest trajectory(start from lowest structure), overall barriers and akmc total time for each transition
"""
#!/usr/bin/env python

import sys
import os
import glob
import pandas as pd
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read, write
from ase.io.trajectory import Trajectory

dynamics = pd.read_table('dynamics.txt', delimiter = r'\s+', skiprows = [0,1], names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])

for i in range(len(dynamics['step-number'])-1):
   if dynamics['product-id'].iloc[i] != dynamics['reactant-id'].iloc[i+1]:
      if isinstance(dynamics['barrier'].iloc[i+1], int): 
         print i, "superbasin"
         continue
      print i, dynamics['barrier'].iloc[i+1]
