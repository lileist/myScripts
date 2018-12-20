#!/usr/bin/env python
"""
This code is used to get Job status, #of fcs and time for each job under debug_results
"""

import sys
import os
import re
import glob
import ast
import pandas as pd
#import matplotlib.pyplot as plt
#from pele.storage import Database
#from pele.landscape import database2graph
#from pele.utils.disconnectivity_graph import DisconnectivityGraph
import numpy
from ase.neighborlist import neighbor_list as nl
from ase.io import read
import collections


arg = sys.argv
current = os.getcwd()
cycle = int(arg[1])
jobid = int(arg[2])
if os.path.isfile('jobstatus'):
   os.system('rm jobstatus')
state_main_dir = current+"/debug_results/"
os.chdir(state_main_dir)

#read in states with energy as dictionary
#with open('state_table') as f:
#   states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])

job_listdir=[]
for f in glob.glob('*'):
   try:
      job_listdir.append(f)
   except:
      continue
job_listdir.sort()

time = 0
succ_job = []
fail = []
in_same_state = []
negative_barrier = []
t_fcs = 0
t_time = 0
t_job = 0
for dir in job_listdir:
   if int(dir.split('_')[0]) > cycle:
      break
   if int(dir.split('_')[1]) > jobid:
      continue
   os.chdir(state_main_dir+str(dir))
   
   t_job += 1
 
   f = open('ll_out', 'r')
   try:
     f_results = open('results.dat','r')
   except:
     continue
   for line in f_results.readlines():
       if 'total_force_calls' in line:
           fcs = int(line.split()[0])
           t_fcs += fcs
       if 'potential_energy_reactant' in line:
           e_rs = float(line.split()[0])
       if 'potential_energy_saddle' in line:
           e_ts = float(line.split()[0])
   for line in f.readlines():
       if 'real' in line:
         time = float(line.split()[1])
         t_time += time
   f.close()

   f = open('ll_out', 'r')
   for line in f.readlines():
       if all( k in line for k in ['Success', 'Final status']):
          succ_job.append([dir, fcs, time, e_ts-e_rs])
          continue
       if  'Saddle is not connected'  in line:
          if 'both minima are the initial state' in line:
             in_same_state.append([dir, fcs, time])
          else:
             fail.append([dir,fcs, time])
          continue
       if any(k in line for k in ['No forward barrier was found', 'Final status: Negative barrier detected']):
          negative_barrier.append([dir,fcs, time])
          continue
       if any(k in line for k in ['Final status: Barrier too high','Final status: Too many iterations','Minimizations from saddle did not converge']):
          fail.append([dir,fcs, time])
          continue
   f.close()


def output_write(output, results, label):
   output.write("%s\n"%(label))
   for item in results:
      try:
        output.write("  %6s %12d %12.2f %12.6f\n"%(item[0], item[1], item[2], item[3]))
      except:
        output.write("  %6s %12d %12.2f\n"%(item[0], item[1], item[2]))
os.chdir(current)
output = open('jobstatus', 'a')
output_write(output, succ_job, 'Success jobs')
output_write(output, fail, 'Failed jobs')
output_write(output, in_same_state, 'In_same_state jobs')
output_write(output, negative_barrier, 'Negative_barrier jobs')

print t_job, t_time, t_fcs
output.write("   succ_ratio  time/success   time/job   fcs/success    fcs/job \n")
try:
  output.write("%8.4f  %14.2f   %14.2f  %14.2f  %14.2f  %14.2f\n"%(float(len(succ_job))/t_job,  float((t_time)/len(succ_job)),float((t_time)/t_job),
                                           float(t_fcs)/len(succ_job), float(t_fcs)/t_job, float(t_time)/t_fcs))
except:
  output.write("%8.4f  %14s   %14.2f  %14s  %14.2f\n"%(float(len(succ_job))/t_job,  "NaN",float((t_time)/t_job),
                                           "NaN", float(t_fcs)/t_job))
