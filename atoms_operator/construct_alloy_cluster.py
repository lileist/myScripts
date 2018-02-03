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
  numberOfElement2 = 163
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
    numb_1 = int(paras['numberOfElement1'])
    numb_2 = int(paras['numberOfElement2'])
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
       indices = random.sample(atom_indices, numb_1)
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
    #define core
    """
      core_define = 1
      max_core = 500    #number of atoms in core
      core_element_1 = 40  #number of element 1 in core
    """
    core_define = int(paras['core_define'])
    if core_define == 1:
       core = []
       core_numb = 0
       core_element_1 = int(paras['core_element_1'])
       max_core = core_element_1 + numb_2
       gr_numb = 0
       for i in range(n):
           if core_numb == max_core:
              break
           gr_numb += len(gr[i])
           for j in range(len(gr[i])):
             core.append(gr[i][j])
             core_numb += 1
             if core_numb == max_core:
                break
       indices = random.sample(core, core_element_1)
       #set core
       for i in core:
           p1[i].symbol = element_2
       print element_1, "in core:",len(indices)
       for i in indices:
           p1[i].symbol = element_1
       #assign others
       shell_1 = 0
       for atom in p1:
           if atom.index not in core:
              p1[atom.index].symbol = element_1
              shell_1 += 1
       print element_1, "in shell:",shell_1
       write(paras['output_structure'], p1,format = paras['output_format'])
       sys.exit()

    for i in range(n):
        if species_numb == numb_1:
           break
        for j in range(len(gr[i])):
            p1[gr[i][j]].symbol = element_1
            species_1.append(gr[i][j])
            species_numb += 1
            if species_numb == numb_1:
               break
    if species_numb < numb_1:
       print "increase n"
    for atom in p1:
        if atom.index not in species_1:
           atom.symbol = element_2
    write(paras['output_structure'], p1,format = paras['output_format'])
if __name__ == '__main__':
    main()
