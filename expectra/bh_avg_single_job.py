#!/usr/bin/env python

import sys
import os
import ase,numpy
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
"""
count_lbfgs [rawNumberOfElement] [MaxBH_Step]
"""
def main():
    directory = os.getcwd()
    arg = sys.argv
    #number = int(arg[2])
    elem_number = int(arg[1])
    output = open("dr.dat",'w')
    total_dr = 0
    total_step = 0
    dr_list = []
    avg_dr = []
    acc_list = []
    pseu_acc_list = []
    force_calls = []
    maxbh_step = int(arg[2])
    force_calls.append(0)
    dr = 0.0
    step = 0
    acc_numb = 0
    pseu_acc = 0
    bh_file='./pot'

    f_bh = open(bh_file+'/pot_log.dat','r')
    f_out = open('out','r')
    bh_lines = f_bh.readlines()
    out_lines = f_out.readlines()

    force_calls = 0
    for line in out_lines:
       if 'iteration' in line:
          force_calls += int(line.split()[1])

    for line in bh_lines:
       if line.startswith("#"):
          continue
       if len(line.split()) < 11:
          break
       fields = line.split()
       curr_pot = float(fields[7])
       u_min = fields[8]
       bh_step = fields[1]
       if fields[1]=='-1':
          pot = curr_pot
       pot_var = curr_pot - pot
       if fields[2]=='1':
          pseu_acc += 1
          pot = curr_pot
          if pot_var != 0.0:
             acc_numb += 1
       dr += float(line.split()[elem_number])
       dr_list.append(float(line.split()[elem_number]))
       step += 1
       if int(bh_step) == maxbh_step:
          break
    fields = bh_lines[len(bh_lines)-1].split()
    pseu_acc_list.append(float(pseu_acc)/float(step))
    acc_list.append(float(acc_numb)/float(step))
    total_dr += dr
    total_step += step
    output.write(" %8.6f %8.6f %s %8.6f %s\n"%(float(pseu_acc)/float(step), float(acc_numb)/float(step), bh_step, force_calls, u_min))
    avg_dr.append(dr/float(step))
    output.write("%s  %8.6f %8.6f\n"%("average", numpy.mean(pseu_acc_list), numpy.mean(acc_list)))
    output.write('%15.6f %15.6f %15.6f %15.6f\n'%(numpy.std(numpy.array(dr_list), ddof=1), numpy.std(numpy.array(dr_list), ddof=0), max(avg_dr), min(avg_dr)))
    output.write(""%())
if __name__ == '__main__':
    main()
    
