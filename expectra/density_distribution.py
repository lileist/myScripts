#!/usr/bin/env python

import sys
import os
#import statistics
import numpy
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
"""
analyze basin_hopping benchmark jobs to figure out statistic data (bh stesps and force calls used, success rate):
count_lbfgs [start_number] [#ofruns] [#OFColumnAccessed] [min] [max] [de] 
"""
def main():
    directory = os.getcwd()
    arg = sys.argv
    number = int(arg[2])
    #column = int(arg[3])
    de = float(arg[3])
    log_e = open('gr.dat', 'w')
    e_gr=[]
    energy=[]
    pot_min = None
    pot_max = None
    for i in range(int(arg[1]), int(arg[1])+number):
       bh_file=directory+"/run-"+str(i)+'/pot'
       #f_bh = open(directory+"/run-"+str(i)+'/pot/pot_log.dat', 'r')
       if os.path.exists(bh_file):
          f_bh = open(bh_file+'/pot_log.dat','r')
       else:
          #num_bh[str(i)] = 0
          print 'run-',i, 'is not run'
          continue
       bh_lines = f_bh.readlines()
       line_number = 0
       for line in bh_lines:
           line_number += 1
           if line.startswith('#'):
              continue
           fields = line.split()
           pot = float(fields[7])
           energy.append(pot)
           if pot_min is None or pot_min > pot:
              pot_min = pot
           if pot_max is None or pot_max < pot:
              pot_max = pot
           if line_number > 50001:
              break
    for i in xrange (int((pot_max - pot_min)/de)+1):
        e_gr.append(0.0)
    for pot in energy:
        for j in xrange (len(e_gr)):
            if pot >= pot_min + float(j)*de and pot < pot_min + float((j+1))*de:
               e_gr[j] += 1
    log_e.write('# %s %s\n'%("std","mean"))
    log_e.write('# %15.6f %15.6f\n'%(numpy.std(numpy.array(energy)), numpy.mean(numpy.array(energy))))
    for i in xrange (len(e_gr)):
        e_gr[i] = float(e_gr[i])/len(energy)
        r = pot_min + float(i)*de - de*0.5
        log_e.write('%15.6f  %15.6f\n' % (r, e_gr[i]))

if __name__ == '__main__':
    main()
    
