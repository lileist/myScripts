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
from math import pi
from copy import deepcopy
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

def read_atoms(filename, state_number = None, mode = None):
    f = open(filename, 'r')
    atoms=[]
    cycle = -1
    if mode == 'wholefile':
       lines = f.readlines()
       count = 0
       count_atom = 0
       elements = []
       positions = []
       for line in lines:
           count += 1
           if count == 1:
              atom_numb = int(line.split()[0])
              continue
           if count <= (atom_numb+2)*state_number+2:
              continue
           fields = line.split()
           elements.append(fields[0])
           positions.append( [ float(fields[j+1]) for j in range(3) ] )
           count_atom += 1
           if count_atom == atom_numb:
              elements = numpy.array(elements)
              positions = numpy.array(positions)
              return Atoms(elements, positions=positions)
    else:
       while True:
           elements = []
           positions = []
           line = f.readline()
           if not line:
              break
           if cycle == -1:
              atom_numb = int(line.split()[0])
           line = f.readline()
           if state_number is not None:
              if cycle == state_number:
                 for i in range (atom_numb):
                     line = f.readline()
                     fields = line.split()
                     elements.append(fields[0])
                     positions.append( [ float(fields[j+1]) for j in range(3) ] )
                 elements = numpy.array(elements)
                 positions = numpy.array(positions)
                 return Atoms(elements, positions=positions)
              else:
                 for i in range (atom_numb):
                     f.readline()
           else:
              for i in range (atom_numb):
                  line = f.readline()
                  fields = line.split()
                  elements.append(fields[0])
                  positions.append( [ float(fields[j+1]) for j in range(3) ] )
              elements = numpy.array(elements)
              positions = numpy.array(positions)
              atoms.append(Atoms(elements, positions=positions))
           cycle += 1
    f.close()
    return atoms

def find_cycle(atoms, nn_waters, waters, end_atom, visited_waters,water_numb=0):
    sys.setrecursionlimit(2000)
    print 'check input', end_atom
    #check if any next nearest water connecting to end atom    
    temp_visited = deepcopy(visited_waters)
    water_numb += 1
    for water in nn_waters:
        if len(temp_visited)==0:
           visited_waters.append(water)
        else:
           for visited_water in temp_visited:
              if water != visited_water:
                 visited_water.append(water)
              else:
                 nn_waters.pop(water)
        for i in range(2):
            if atoms.get_distance(end_atom, water[i+1].index) <1.8 and atoms.get_angle([end_atom, water[0].index, water[i+1].index])*180/pi < 30:
               return water_numb
    #if none conneting to end atom, find next nearest water. The current water itself needs to be ruled out
    for nn_water in nn_waters:
        nnn_waters=[]
        for water in waters:
            #rule out current water and those visisted in previous cycle
            if water == nn_water or water in visited_waters or water in nnn_waters:
               continue
            for i in range(2):
                print nn_water[i+1].index, nn_water[0].index, water[0].index
                if atoms.get_distance(nn_water[i+1].index, water[0].index) <1.8 and atoms.get_angle([nn_water[i+1].index, nn_water[0].index, water[0].index])*180/pi < 30:
                   nnn_waters.append(water)
                   break
                if atoms.get_distance(nn_water[0].index, water[i+1].index) <1.8 and atoms.get_angle([water[i+1].index, water[0].index, nn_water[0].index])*180/pi < 30:
                   nnn_waters.append(water)
                   break
    #all waters have been considered
    print 'nnn:',len(nnn_waters)
    if len(nnn_waters)==0:
       return water_numb
    find_cycle(atoms, nn_waters=deepcopy(nnn_waters), waters=deepcopy(waters), end_atom=end_atom, visited_waters=deepcopy(visited_waters),water_numb=water_numb)
                  
        

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    configs = read_atoms(paras['configs_file'])
    start_atom= [int(field) for field in paras['start_atom'].split()]
    end_atom = int(paras['end_atom'])
    non_list = [int(field) for field in paras['non_list'].split()]
    image_numb = 0
    output = open(paras['outputfile'], 'w')
    #group atom: H2O, NH3
    water_numb = 0
    for atoms in configs:
        image_numb+=1
        waters = []
        for atom in atoms:
            if atom.symbol == 'O' and atom.index not in non_list:
               waters.append([atom,0,0,0])
            if atom.symbol == 'H':
               for water in waters:
                   if atoms.get_distance(water[0].index, atom.index) < 1.2:
                      if water[3] == 0:
                         water[1] = atom
                         water[3] += 1
                         continue
                      water[2] = atom
                      water[3] += 1
        print waters
        for i in start_atom:
            for water in waters:
                if atoms.get_distance(i, water[0].index) < 2.9:
                   water_numb=find_cycle(atoms, nn_waters=[water], waters=deepcopy(waters), end_atom=end_atom, visited_waters=[], water_numb=0)
        output.write('%d %d'%(image_numb, water_numb))
if __name__ == '__main__':
    main()
    
