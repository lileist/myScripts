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
    
def find_distance(filename, index_1, index_2):
    f = open(filename,'r')
    line_numb = 0
    distance = []
    while True:
       line = f.readline()
       if not line:
          break
       line_numb+=1
       if line_numb == 1:
          atom_numb = int(line.split()[0])
       line = f.readline()
       positions = []
       for i in range(atom_numb):
           line = f.readline()
           if i == index_1 or i == index_2:
              fields = line.split()
              positions.append( numpy.array([ float(fields[j+1]) for j in range(3) ] ))
       #distance.append(la.norm(positions[0]-positions[1]))
       distance.append(math.sqrt((positions[0][0]-positions[1][0])**2+(positions[0][1]-positions[1][1])**2+(positions[0][2]-positions[1][2])**2))
    f.close()
    return numpy.array(distance)

def main():
    arg = sys.argv
    input = open(arg[1], 'r')
    output = open('countour.dat','w')
    x_numb = 0
    while True:
        for i in range(100):
            line = input.readline()
            if not line:
               break
            fields = line.split()
            output.write("%s "%(fields[2]))
        x_numb+=1
        output.write("\n")
        if not line:
           break
    print x_numb
if __name__ == '__main__':
    main()
    
