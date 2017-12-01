#!/usr/bin/env python

import sys
import os
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
"""
analyze basin_hopping benchmark jobs to figure out statistic data (bh stesps and force calls used, success rate):
count_lbfgs [start_number] [#ofruns] [stopCriterion] [#ofColumnValuesAccesed] [#maxBHsteps]
"""
def main():
    directory = os.getcwd()
    arg = sys.argv
    number = int(arg[2])
    column_numb = int(arg[4])
    num_bh = {}
    num_fc = {}
    converged = 0
    incompleted = []
    total_bh = 0
    total_fc = 0
    print float(arg[3])
    if len(arg)<6:
       max_bhsteps=50000
    else:
       max_bhsteps=int(arg[5])
    for i in range(int(arg[1]), int(arg[1])+number):
       bh_file=directory+"/run-"+str(i)+'/pot'
       #f_bh = open(directory+"/run-"+str(i)+'/pot/pot_log.dat', 'r')
       geo_opt = directory+"/run-"+str(i)+'/geo_opt.log'
       if os.path.exists(bh_file):
          f_bh = open(bh_file+'/pot_log.dat','r')
       else:
          #num_bh[str(i)] = 0
          print 'run-',i, 'is not run'
          continue
       if os.path.exists(geo_opt):
          f_fc = open(geo_opt, 'r')
          fc_lines = f_fc.readlines()
       else:
          fc_lines = []
       bh_lines = f_bh.readlines()
       num_bh[str(i)] = (sum(1 for line in bh_lines))
       num_fc[str(i)] = (sum(1 for line in fc_lines))
       if num_fc[str(i)] == 1:
          num_fc[str(i)] = int(fc_lines[0].split()[0])
       fields = bh_lines[len(bh_lines)-1].split()
       print len(bh_lines)
       if float(fields[column_numb]) <= float(arg[3]) and len(bh_lines) < max_bhsteps+1:
          converged+=1
          total_bh += num_bh[str(i)]
          total_fc += num_fc[str(i)]
       else:
          print float(fields[column_numb])-float(arg[3])
          num_bh.pop(str(i))
          num_fc.pop(str(i))
          if int(fields[1]) <= int(max_bhsteps-1):
            incompleted.append([i,fields[1]])
       f_bh.close()
       if os.path.exists(geo_opt):
          f_fc.close()
    if converged == 0:
       print "no job coverged"
    else:
       print "bh steps:",float(total_bh)/float(converged), num_bh[max(num_bh, key=lambda i:num_bh[i])], num_bh[min(num_bh, key=lambda i:num_bh[i])]
       print "fc steps:",float(total_fc)/float(converged), num_fc[max(num_fc, key=lambda i:num_fc[i])], num_fc[min(num_fc, key=lambda i:num_fc[i])]
    print "not finised jobs:", incompleted
    print "success rate:",float(converged)/float(number)

if __name__ == '__main__':
    main()
    
