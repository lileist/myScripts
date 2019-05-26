#!/usr/bin/env python
"""
This code is used to integrate the mean force from thermal integration, thus obtaining free energy
datAlz_TherInteg.py [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    cp2k_inp = suppl.inp
    job_submit_script = qsub.hf
    job_submit_cmd    = sbatch
"""
import sys
import os
import errno
import ase
from ase import Atoms
from ase.io import read, write, Trajectory
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
#from ase.constraints import constrained_indices, FixAtoms
from numpy import linalg as la
import numpy, math
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

    
def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    distances=paras['distance'].split()
    md_traj = read('3w-pos-1.xyz',index=slice(-3000,None),format='xyz')
    for p in md_traj:
       dist = p.get_distance(7,13) + p.get_distance(9,11)
       for distance in distances:
           if dist > float(distance)-0.05 and dist < float(distance)+0.05:
              make_dir(distance)
              write(distance+'/geometry.xyz', p, format='xyz')
              break
    traj = Trajectory('selected.traj','w')
    for distance in distances:
        traj.write(read(distance+'/geometry.xyz', index=':', format='xyz')[0])
if __name__ == '__main__':
    main()
    
