#!/usr/bin/env python
"""
This code is used to integrate the mean force from thermal integration, thus obtaining free energy
datAlz_TherInteg.py [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
distance = 2.4  2.5  2.6  2.8  3.0  3.2  3.4 3.5 3.6  3.8  4.0  4.2  4.4  4.6  4.8 5.0
mode = accumulative or separate
cp2k_inp = suppl.inp
job_submit_script = run-cp2k.sub
job_submit_cmd    = sbatch
stablize_step = 20000
d_step = 4000
max_step = 400000
reactant = 4.2  5.0
product  = 2.4  2.8
dx = 0.1
index_1 = 14
index_2 = 16
"""
import sys
import os
import errno
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
#from ase.constraints import constrained_indices, FixAtoms
from collections import OrderedDict 
from numpy import linalg as la
import numpy as np
import math
import subprocess

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
           raise

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

def format_joint(elements):
    jointed = ''
    for element in elements:
       jointed+="{:12.6f}".format(element)
    return jointed


def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    dirs=paras['dirs'].split()
    target_file = paras['target_file']

    output = open(dirs[0].split('_')[0]+'_'+target_file.split('_')[0][0:4]+'_std.dat','w')

    output.write('{:12s} {:12s} {:12s} {:12s} {:12s} {:12s} {:12s}\n'.format(
                   'T', 'reactionE','avg','std', 'barrier', 'std', 'avg'))
    ess = {}
    for d in dirs:
       ess_input = open('../'+d+'/'+target_file, 'r')
       for line in ess_input.readlines():
         if 'time' in line:
           ess[d] = []
           continue
         ess[d].append([float(field) for field in line.split()])
       ess_input.close()
       stds = np.std(ess[d], axis=0, ddof=1)
       avgs = np.average(ess[d], axis=0)
       output.write('{:12s} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n'.format(
                    d.split('_')[1][0:-1], ess[d][-1][4], avgs[4], stds[4], ess[d][-1][5], avgs[5], stds[5]))

if __name__ == '__main__':
    main()
    
