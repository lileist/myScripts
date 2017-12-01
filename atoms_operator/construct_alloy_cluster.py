#!/usr/bin/env python
import sys
import os
import random
import shutil
from ase.io import Trajectory
from ase.io import read,write
from ase import Atoms
from numpy import linalg as la
""" 
constructure alloy clusters (core shell or random):
input example:
  input_structure = POSCAR
  input_format = vasp
  output_structure = CONTCAR
  output_format = vasp
  random = 1                 #generate structure randomly
  center = cm
  numberOfElement1 = 163
  rmin = 1.0
  dr = 0.5
  numberOfInterval = 6
  element_1 = 'Pt'
  element_2 = 'Au'
"""
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
    p1 = read(paras['input_structure'], index=0, format=paras['input_format'])
    if paras['center']== 'cm':
       center = p1.get_center_of_mass()
    else:
       center= [float(field) for field in paras['center'].split()]
    element_1 = str(paras['element_1'])
    element_2 = str(paras['element_2'])
    numb = int(paras['numberOfElement1'])
    rmin = float(paras['rmin'])
    dr = float(paras['dr'])
    n = int(paras['numberOfInterval'])
    random_stru = int(paras['random'])
    species_numb = 0
    species_1 = []
    numb_reached = False
    gr = []
    #generate structure randomly
    if random_stru == 1:
       atom_indices = range(0, len(p1))
       indices = random.sample(atom_indices, numb)
       for atom in p1:
           assigned = False
           for index in indices:
               if atom.index == index:
                  atom.symbol = element_1
                  assigned = True
           if not assigned:
                  atom.symbol = element_2
       write(paras['output_structure'], p1,format = paras['output_format'])
       sys.exit()

    for i in range(n):
        gr.append([])
    #analyze radical distribution
    for atom in p1:
        if numb_reached:
           break
        vect = atom.position - center
        dist = la.norm(vect)
        for i in range(n):
            if dist >= rmin + float(i)*dr and dist < rmin +float(i+1)*dr:
               gr[i].append(atom.index)

    for i in range(n):
        if species_numb == numb:
           break
        for j in range(len(gr[i])):
            p1[gr[i][j]].symbol = element_1
            species_1.append(gr[i][j])
            species_numb += 1
            if species_numb == numb:
               break
    if species_numb < numb:
       print "increase n"
    for atom in p1:
        if atom.index not in species_1:
           atom.symbol = element_2
    write(paras['output_structure'], p1,format = paras['output_format'])
if __name__ == '__main__':
    main()
