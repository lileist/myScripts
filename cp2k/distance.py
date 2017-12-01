#!/usr/bin/env python
"""
This code is used to calculate strain of nanoparticles.
Method: stretch partilces along certain directions by dx and fix the boundary atoms, do optimization
strain_cal_cluster [inputfile]
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
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
from ase.constraints import constrained_indices, FixAtoms
import numpy
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
    cp2k_inp = open(paras['cp2k_inp'], 'r')
    lines = cp2k_inp.readlines()
    distances=paras['distance'].split()
    print distances
    for distance in distances:
        make_dir(distance)   
        output = open(distance+'/suppl.inp', 'w')
        for line in lines:
            if all(s in line for s in ["TARGET",  "[angstrom]"]):
               fields = line.split('[angstrom]')
               output.write("%s %s %s\n"%(fields[0], '[angstrom]', distance))
               continue
           # Used to remove all empty lines
           # if len(line.split())<1:
           #    continue
            output.write("%s"%(line))
        os.system('cp geometry.xyz '+distance)
        os.system('cp '+paras['job_submit_script']+' '+distance)
        with cd(distance):
             os.system(paras['job_submit_cmd']+' '+paras['job_submit_script'])
             print os.getcwd(),'is submitted'
        print 'Back to directory',os.getcwd()
if __name__ == '__main__':
    main()
    
