#!/usr/bin/env python
"""
Code used to analyze DOSCAR and output pdos for the selected atoms
pdos.py [atoms_selected]
Files required: POSCAR DOSCAR PROCAR
Two ways to define the selected atoms:
    1. provide the index list of the atoms you're interested:
       pdos.py 121 122
    2. provide the atom type and a box which is defined by minimum and maximum xyz.
       pdos.py Ti 
       The code will ask you to input minimum and maximum xyz values to define the box
"""
from ase import Atoms
from ase.io import read
import sys
import os
import numpy as np
import copy

def pick_atoms(atoms, atom_type, 
               x_max = 100.0, 
               x_min = -100.0, 
               y_max = 100.0, 
               y_min = -100.0, 
               z_max = 100.0, 
               z_min = -100.0):
    selected_atoms = []
    for atom in atoms:
       if atom.symbol in atom_type:
          if atom.x < x_max and atom.x > x_min:
             if atom.y < y_max and atom.y > y_min:
                if atom.z < z_max and atom.z > z_min:
                   #print atom.index
                   selected_atoms.append(atom.index)
    return selected_atoms

def readBaderCharge(filename)
    f = open(filename, 'r')
    while True:
          line = f.readline()
          if not line:
             break

          
def assign_min_max(inputs):
    fields = [ float(field) for field in inputs.split()]
    if len(fields) > 0:
       min = fields[0]
       max = fields[1]
    else:
       min = -1000.0
       max = 1000.0
    return min, max
           
def main():
    arg = sys.argv
    atom_type = []
    for i in range (len(arg)):
        if i > 0 :
           atom_type.append(arg[i]) 
    index = True
    try:
       #atoms selected based on provided index
       selected_atoms = [int(field) for field in atom_type]
    except:
       index = False
       pass
    #atoms selected based on defined box and atom_type
    if not index:
       x = raw_input("please enter x min and max\n")
       y = raw_input("please enter y min and max\n")
       z = raw_input("please enter z min and max\n")
       x_min, x_max = assign_min_max(x)
       y_min, y_max = assign_min_max(y)
       z_min, z_max = assign_min_max(z)
       print x_min, x_max, y_min, y_max, z_min, z_max
       p1 = read(filename = 'POSCAR', index=0, format = 'vasp')
       selected_atoms = pick_atoms(p1, atom_type, x_max, x_min, y_max, y_min, z_max, z_min)
    print selected_atoms

    pdosfile = open( ''.join(atom_type) + '_pdos.dat', 'w')
    print pdosfile
    if sum(d_dos) == 0.0:
       d_center = 0.0
    else:
       d_center = np.dot(d_dos, np.array(energy)) / sum(d_dos)
    pdosfile.write("%s %8.3f %8.3f %8.3f %8.3f\n" % ("#", fermi, s_center, p_center, d_center))

    pdosfile.close()
if __name__ == '__main__':
    main()
