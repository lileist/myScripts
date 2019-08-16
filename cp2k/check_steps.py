#!/usr/bin/env python
"""
This code is used to integrate the mean force from thermal integration, thus obtaining free energy
datAlz_TherInteg.py [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
cp2k_inp = suppl.inp
job_submit_script = run-cp2k.sub
job_submit_cmd    = sbatch
distance = 2.4  2.5  2.6  2.8  3.0  3.2  3.4 3.5 3.6  3.8  4.0  4.2  4.4  4.6  4.8 5.0
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

def simpson(f, a, b, n):
    if n%2:
       raise ValueError("n must be even (received n=%d)" % n)
    h = (b - a) / n
    s= f[str(a)] +f[str(b)]
    
    for i in range(1, n, 2): 
        s += 4 * f[str(a + i * h)]
    for i in range(2, n-1, 2):
        s += 2 * f[str(a + i * h)]
    return s * h / 3
    
def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    d_step = int(paras['d_step'])
    stablize_steps = int(paras['stablize_step'])
    distances=paras['distance'].split()
    
    rs = paras['reactant'].split()
    fs = paras['product'].split()
    ts = paras['transition_state'].split()
    rs_min = distances.index(rs[0])
    rs_max = distances.index(rs[1])
    fs_min = distances.index(fs[0])
    fs_max = distances.index(fs[1])
    ts_min = distances.index(ts[0])
    ts_max = distances.index(ts[1])


    dx = float(paras['dx'])
    forces=OrderedDict()
    avg_forces = OrderedDict()
    free_energies = OrderedDict()
    md_steps = []
    atom_distance = {}
    force_file = None
    for distance in distances:
        forces[distance]=[]
        if force_file is None:
           cp2k_inp = open(distance+'/'+paras['cp2k_inp'], 'r')
           lines = cp2k_inp.readlines()
           for line in lines:
               if 'PROJECT_NAME' in line:
                  force_file = line.split()[1]+'-1.LagrangeMultLog'
                  break
        input_force = open(distance+'/'+force_file, 'r')
        total_force = 0
        shake_numb = 0
        numb_force = 0
        numb_line = 0
        lines = input_force.readlines()
        for line in lines:
#            numb_line += 1
            if 'Shake  Lagrangian' in line:
               forces[distance].append(float(line.split()[3]))
        forces[distance]=np.array(forces[distance])
        input_force.close()
        n_forces = len(forces[distance])
        md_steps.append(n_forces)
    print md_steps
    n_slices = int((np.amin(md_steps)-stablize_steps)/d_step)+1
    print n_slices
if __name__ == '__main__':
    main()
    
