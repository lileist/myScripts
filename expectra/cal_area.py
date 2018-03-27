#!/usr/bin/env python
import numpy as np
import sys

def calc_area(y_exp, y_theory, calc_type='area', average = False):
    if len(y_exp) != len(y_theory):
        print "Warning: number of points in chi_exp and chi_theory is not equal"
        
    numb = min(len(y_exp),len(y_theory))
    area_diff = 0.00
    for i in range(0, numb):
      diff = np.absolute(y_exp[i] - y_theory[i])
      if calc_type == 'least_square':
        area_diff = area_diff + (diff/y_exp[i])**2
      elif calc_type == 'area':
        area_diff = area_diff + diff
    if average:
       area_diff = area_diff / numb
    #print ('%s: %15.6f' % ("area_diff", area_diff))
    return area_diff

arg = sys.argv
y1 =np.loadtxt(arg[1])
y2 =np.loadtxt(arg[2])
print calc_area(y1[:,1], y2[:,1])

