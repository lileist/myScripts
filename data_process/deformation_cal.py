#!/usr/bin/env python
"""
This code is used to compare two structures and calculate the change of each bond and output in a matrix format
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
import numpy as np

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

def cal_distance(atoms, atom_index):
    distance_lib = []
    for i in range(len(atom_index)):
        if i == len(atom_index)-1:
           break
        for j in range(i+1, len(atom_index)):
            ij_dist = atoms.get_distance(atom_index[i], atom_index[j])
            if ij_dist < 3.4:
              position = (atoms[atom_index[i]].position + atoms[atom_index[j]].position) / 2
              distance_lib.append([position, ij_dist, np.array([atom_index[i],atom_index[j]])])
    return np.array(distance_lib)

def printMatrix(a, e_min, e_max, s_min, s_max):
   print "# Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]"
   print "# "+ ("%15.6f" % e_min) +("  %15.6f" % e_max) +("  %15.6f" % s_min) +("  %15.6f" % s_max) 
   rows = a.shape[0]
   cols = a.shape[1]
   for i in range(0,rows):
      for j in range(0,cols):
         print("%8.5f" %a[i,j]),
      print
   print  

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    p1 = read(paras['configs_1'], index=0, format='vasp')
    p2 = read(paras['configs_2'], index=0, format='vasp')
    t_deform = 0
    #find index of target species
    atom_index = []
    p3 = Atoms()
    for i in range(len(p1)):
        if p1[i].symbol in paras['target_species']:
           atom_index.append(i)
           p3.append(p1[i])
    #print atom_index
    min_x = np.min(p3.get_positions()[:,0])
    max_x = np.max(p3.get_positions()[:,0])
    min_y = np.min(p3.get_positions()[:,1])
    max_y = np.max(p3.get_positions()[:,1])
    #create grids
    dr = 0.5
    grids = np.zeros((int((max_y-min_y)/dr)+1, int((max_x-min_x)/dr)+1))
    distances_1 = cal_distance(p1, np.array(atom_index))
    #distances_2 = cal_distance(p2, np.array(atom_index))
    #print distances_1
    #print distances_2.shape
    dist_max = 0
    bond_n=0
    for i in range(len(distances_1)):
       #print distances_1[i][2]
       new_dist = p2.get_distance(distances_1[i][2][0],distances_1[i][2][1])
       bond_n += 1
       #distances_1[i][1]=distances_2[i][1]-distances_1[i][1]
       for k in range(int((max_y-min_y)/dr)+1):
          for j in range(int((max_x-min_x)/dr)+1):
            if distances_1[i][0][0] >= min_x + j*dr and distances_1[i][0][0] < min_x + (j+1)*dr \
               and distances_1[i][0][1] >= min_y + k*dr and distances_1[i][0][1] < min_y + (k+1)*dr:
               grids[k][j] = new_dist-distances_1[i][1]
               t_deform += grids[k][j]
               if grids[k][j] > dist_max:
                  dist_max = grids[k][j]
                  largest=distances_1[i][2]

    print largest, dist_max
    print t_deform, t_deform/bond_n
    printMatrix(np.array(grids), min_x, max_x, min_y, max_y)

if __name__ == '__main__':
    main()
    
