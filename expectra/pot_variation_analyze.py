#!/usr/bin/env python

import sys
import os
import ase,math
import numpy
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
"""
analyze basin_hopping benchmark jobs to check variation of assessed value. gr function and difference of the value relative to current accepted value will be recorded :
count_lbfgs [start_number] [#ofruns] [de]
"""
def main():
    directory = os.getcwd()
    arg = sys.argv
    number = int(arg[2])
    #column_numb = int(arg[4])
    de = float(arg[3])
    pot_var_list = []
    pot_min = None
    pot_max = None
    for i in range(int(arg[1]), int(arg[1])+number):
       output = open('pot_var'+str(i)+'.dat','w')
       bh_file=directory+"/run-"+str(i)+'/pot'
       if os.path.exists(bh_file):
          f_bh = open(bh_file+'/pot_log.dat','r')
       else:
          #num_bh[str(i)] = 0
          print 'run-',i, 'is not run'
          continue
       bh_lines = f_bh.readlines()
       #num_bh[str(i)] = (sum(1 for line in bh_lines))
       line_number = 0
       for line in bh_lines:
          line_number += 1
          if line.startswith('#'):
             continue
          fields = line.split()
          curr_pot = float(fields[7])
          if fields[1]=='-1':
             pot = curr_pot
          pot_var = curr_pot - pot
          if fields[2]=='1':
             pot = curr_pot
          output.write("%s  %15.6f\n"%(fields[1],pot_var))
          pot_var_list.append(pot_var)
          if fields[1] == '-1':
             continue
          if pot_min is None or pot_min > pot_var:
             pot_min = pot_var
          if pot_max is None or pot_max < pot_var:
             pot_max = pot_var
          if line_number > 50001:
             break
    #get gr fuction
    pot_number = 0
    pot_gr = []
    pot_log = []
    pot_positive =[]
    pot_negative=[]
    for pot_var in pot_var_list:
        pot_number += 1
        if pot_var <=0:
            pot_log.append(0.0)
            pot_negative.append(pot_var)
        else:
            pot_log.append(math.log(pot_var))
            pot_positive.append(pot_var)
        for i in xrange (int((pot_max - pot_min)/de) +1 ):
            #initialize pot_gr
            if pot_number == 1:
               pot_gr.append(0)
            if pot_var >= pot_min + float(i) * de and pot_var < pot_min + float(i+1) * de:
               pot_gr[i] += 1
    log_gr = open('pot_var_gr.dat', 'w')
    log_gr.write("#%15.6f %15.6f %15.6f\n"%(pot_max-pot_min, numpy.std(numpy.array(pot_var_list)), numpy.mean(pot_var_list)))
    log_gr.write("#%15.6f %15.6f\n"%(numpy.std(numpy.array(pot_log)), numpy.mean(pot_log)))
    log_gr.write("#%s %s %s %s %s %s %s\n"%("neg_std","neg_mean","pos_std","pos_mean","neg_cv","pos_cv","ratio"))
    neg_std = numpy.std(numpy.array(pot_negative))
    neg_mean=numpy.mean(pot_negative)
    pos_std =numpy.std(numpy.array(pot_positive))
    pos_mean=numpy.mean(pot_positive)
    neg_cv = abs(neg_std/neg_mean)
    pos_cv = pos_std/pos_mean
    log_gr.write("#%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n"%(neg_std, neg_mean, pos_std, pos_mean, neg_cv, pos_cv, neg_cv/pos_cv))
    negative_area = 0.0
    positive_area = 0.0
    for i in xrange (len(pot_gr)):
        pot_gr[i] = float(pot_gr[i])/float(len(pot_var_list))
        r = pot_min + float(i) * de - de*0.5
        log_gr.write('%15.6f %15.6f\n' % (r, pot_gr[i]))
        if r <=0:
           negative_area = negative_area + de*pot_gr[i]
        else:
           positive_area += de*pot_gr[i]
    print "negative area, positive area, ratio:", negative_area, positive_area, negative_area/positive_area
if __name__ == '__main__':
    main()
    
