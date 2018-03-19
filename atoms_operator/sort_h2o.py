#!/usr/bin/env python
"""
This code is used to calculate strain of nanoparticles.
Method: stretch partilces along certain directions by dx and fix the boundary atoms, do optimization
strain_cal_cluster [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    configs_file = 3w-pos-1.xyz
    central_index = 1 2 3
    target_index = 5 6
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

def main():
    arg = sys.argv
    configs = read(arg[1],index=0,format='xyz')
    configs.set_cell([[19.9683,0,0],[0,19.9683,0],[0,0,35.0]],scale_atoms=False,fix=None)
    configs.set_pbc((True, True, True))
    wat=[]
    oxygen=[]
    hydrogen=[]
    nitrogen=[]
    wat_numb = 0
    for atom in configs:
        if atom.symbol == 'O':
          oxygen.append(atom)
        if atom.symbol == 'H':
          hydrogen.append(atom)
        if atom.symbol == 'N':
          nitrogen.append(atom)
    output = open('sorted_water.xyz', 'w')
    for o in oxygen:
        wat_temp = []
        for h in hydrogen:
           if configs.get_distance(o.index,h.index,mic=True) < 1.2:
              if len(wat_temp) == 0:
                 wat_temp.append(o)
              wat_temp.append(h)
        if len(wat_temp)>2:
           wat_numb += 1
           output.write('%s %15.6f %15.6f  %15.6f\n' % (wat_temp[0].symbol, wat_temp[0].x, wat_temp[0].y, wat_temp[0].z))
           output.write('%s %15.6f %15.6f  %15.6f\n' % (wat_temp[1].symbol, wat_temp[1].x, wat_temp[1].y, wat_temp[1].z))
           output.write('%s %15.6f %15.6f  %15.6f\n' % (wat_temp[2].symbol, wat_temp[2].x, wat_temp[2].y, wat_temp[2].z))
    print wat_numb
if __name__ == '__main__':
    main()
    
