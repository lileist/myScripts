#!/usr/bin/env python

import sys
import numpy as np

args = sys.argv
f = open(args[1], 'r')
strain =[]
energy =[]
for line in f.readlines():
   fields = line.split()
   strain.append(float(fields[0]))
   energy.append(float(fields[1]))
curvature_s = np.polyfit(strain, energy, 2)
print curvature_s
