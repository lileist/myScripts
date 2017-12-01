#!/usr/bin/env python
"""
This code is used to analyze the distribution of energy and s (chi_deviation).
How to use:
    gr_ES [inputfile_name] [de] [ds]
"""
import sys
import os
#from expectra.io import read_dots
from expectra.io import read_atoms
import numpy as np

def read_charge(filename, atom_numb, first_nb, second_nb, third_nb):
    f = open(filename)
    first_charge = []
    second_charge = []
    third_charge = []
    while True:
        line = f.readline()
        if not line:
           break
        for i in range(4):
            line = f.readline()
        first = 0.0
        second = 0.0
        third = 0.0
        for i in range(atom_numb):
            line = f.readline()
            fields = [ field for field in line.split()]
            if i in first_nb:
               first = first + float(fields[7])
            if i in second_nb:
               second = second + float(fields[7])
            if i in third_nb:
               third = third + float(fields[7])
        first_charge.append(first)
        second_charge.append(second)
        third_charge.append(third)
        line = f.readline()
        line = f.readline()
        line = f.readline()
    f.close()
    return first_charge, second_charge, third_charge

def log_atoms(first_chg, second_chg, third_chg):
    log_trajectory = open('charge.dat', 'w')
    for i in range(len(first_chg)):
        log_trajectory.write("%d  %15.6f  %15.6f  %15.6f\n" % (i,
                                    first_chg[i], second_chg[i], third_chg[i]))
    log_trajectory.close()

def main():
    arg = sys.argv
    traj_file = arg[1]
    charge_file = arg[2]
    
    cutoff = [3.1, 5.2, 6.5]
    first_nb = []
    second_nb = []
    third_nb = []
    traj = read_atoms(traj_file)
    atoms = traj[0]
    center = [6.36, 8.37, 14.13]
    for atom in atoms:
        if atom.symbol == 'Ce':
           dist = np.linalg.norm(np.array(atom.position-center))
           if dist <= center[0]:
              first_nb.append(atom.index)
           if dist <= center[1] and dist >= center[0]:
              second_nb.append(atom.index)
           if dist <= center[2] and dist >= center[1]:
              third_nb.append(atom.index)
    first_chg, second_chg, third_chg = read_charge(charge_file, len(atoms), first_nb, second_nb, third_nb)   
    log_atoms(first_chg, second_chg, third_chg)

if __name__ == '__main__':
    main()

        
