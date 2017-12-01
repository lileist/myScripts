#!/usr/bin/env python
"""
This code is used to calculate strain of nanoparticles.
Method: stretch partilces along certain directions by dx and fix the boundary atoms, do optimization
strain_cal_cluster [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    configs_file = POSCAR
    target_species = Au
    minimum_distance = 1.0
    gr_number = 100
    outputfile = dist_gr.dat
    step_size = 0.02
"""
import sys
import os
import errno
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
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
    atoms = read(paras['configs_file'], index=0, format='vasp')
    dist_gr = []
    numb = 0
    dist_min = float(paras['minimum_distance'])
    gr_numb = int(paras['gr_number'])
    de = float(paras['step_size'])
    image_numb = 0
    for i in range(gr_numb):
        dist_gr.append(0.0)
    #for atoms in configs:
    atom_index = []
    for i in range(len(atoms)):
        if atoms[i].symbol==paras['target_species']:
           atom_index.append(i)
    print atom_index
    for i in range(len(atom_index)):
        if i == len(atom_index)-1:
           break
        for j in range(i+1, len(atom_index)):
            ij_dist = atoms.get_distance(atom_index[i], atom_index[j])
            for k in range (gr_numb):
                if ij_dist >= dist_min + float(k)*de and ij_dist < dist_min + float((k+1))*de:
                   dist_gr[k] += 1
                   numb+=1
                   break
    output = open(paras['outputfile'], 'w')
    for i in range(gr_numb):
        dist_gr[i] = float(dist_gr[i])/numb
        r = dist_min + float(i)*de - de*0.5
        output.write('%15.6f  %15.6f\n' % (r, dist_gr[i]))
if __name__ == '__main__':
    main()
    
