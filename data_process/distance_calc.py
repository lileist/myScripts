#!/usr/bin/env python

import sys
import os
import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory

"""
Calculate distance between two atoms by reading xyz file
"""
def main():
    arg = sys.argv
    inputfile = arg[1]
    outputfile = arg[2]
    atom1 = int(arg[3])
    atom2 = int(arg[4])
    f = open(inputfile)
    output = open(outputfile,'w')
    atoms = None
    atom_output = Trajectory('trajectory.traj','w',atoms)
    iteration = 0
    trajectories = []
    a = 10.0
    b = 10.0 
    c = 10.0
    cell_atoms = np.array([[a,0,0],[0,b,0],[0,0,c]])
    while True:
       coordinates=[]
       type_atoms = []
       line = f.readline()
       if not line:
          break
       atom_numb = int(line.split()[0])    
       line = f.readline()
       for i in range(0, atom_numb):
           line=f.readline()
           fields = line.split()
           type_atoms.append(fields[0])
           coordinates.append([float(fields[1]),float(fields[2]),float(fields[3])])
       traj = Atoms(np.array(type_atoms), positions=np.array(coordinates), cell=cell_atoms, pbc=True)
       traj.wrap()
       atom_output.write(traj)
       output.write("%6d %8.4f\n"%(iteration, traj.get_distance(atom1, atom2,mic=True)))
       iteration += 1
if __name__ == '__main__':
    main()
    
