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
from ase.io import read, write, Trajectory
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

def main():
    arg = sys.argv
    try:
       configs = read(arg[1], index=":")
    except:
       configs = read_atoms(arg[1])

    traj = Trajectory('selected_'+arg[1],'w')
    chem_symbols = ['H','Pd']
    H_indices  = [atom.index for atom in configs[0] if atom.symbol == 'H']
    Pd_indices = [atom.index for atom in configs[0] if atom.symbol == 'Pd']
    
    #temp = copy.deepcopy(configs)
    for itrj in range(len(configs)):
        print itrj
        config = configs[itrj]
        delete = False
        for i in range(len(H_indices)-1):
           for j in range(i+1, len(H_indices)):
              if config.get_distance(H_indices[i], H_indices[j]) < 0.5:
                 delete = True
                 break
           if delete:
              break

           for pd in Pd_indices:
              if config.get_distance(H_indices[i], pd) < 1.00:
                 delete = True
                 break
           if delete:
              break

        if delete:
           continue
        for i in range(len(Pd_indices)-1):
           for j in range(i+1, len(Pd_indices)):
              if config.get_distance(Pd_indices[i], Pd_indices[j]) < 2.00:
                 delete = True
                 break
           if delete:
              break
        if delete:
           continue
        traj.write(config)

#   deletes.sort(reverse=True)
#   for index in deletes:
#      configs.pop(i=index)
   
   

if __name__ == '__main__':
    main()
    
