#!/usr/bin/env python
"""
This code is used to move target molecules close to target molecule
An example of inputfile (any line starts with '#' will be ignored):
    configs_file = 3w-pos-1.xyz
    distances = 2.2 2.4
    target =  24 26
    destination =  56 58
    geo_center_1 = [int(field) for field in paras['geo_center_1'].split()]
    geo_center_2 = [int(field) for field in paras['geo_center_2'].split()]
"""
import sys
import os
import errno
import ase
import math
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

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    workdir = os.getcwd()
    distances = [str(field) for field in paras['distances'].split()]
    target = [int(field) for field in paras['target'].split()]
    destination = [int(field) for field in paras['destination'].split()]
    geo_center_1 = [int(field) for field in paras['geo_center_1'].split()]
    geo_center_2 = [int(field) for field in paras['geo_center_2'].split()]
    bond_length = 2.0
    output = open('distance.ini','w')
    numb = 0
    for distance in distances:
        os.chdir(workdir+'/../thermal_integ/'+distance)
        os.system('tail -n 807 '+paras['configs_file']+' >'+'temp.xyz')
        atoms = read_atoms('temp.xyz')[0]
        #group atom: H2O, NH3
        v1 = atoms[destination[0]].position - atoms[target[0]].position
        norm_v1 = numpy.linalg.norm(v1)
        atoms[target[0]].position = atoms[target[0]].position + v1* (1-bond_length/norm_v1)
        atoms[target[1]].position = atoms[target[1]].position + v1* (1-bond_length/norm_v1)
        write('geometry_1.xyz',images=atoms)
        geo_1 = numpy.array([0.,0.,0.])
        geo_2 = numpy.array([0.,0.,0.])
        for atom in atoms:
            if atom.index in geo_center_1:
               geo_1 += atom.position
            if atom.index in geo_center_2:
               geo_2 += atom.position
        geo_1 = geo_1 / float(len(geo_center_1))
        geo_2 = geo_2 / float(len(geo_center_2))
        new_distance = round(numpy.linalg.norm(geo_1 - geo_2), 1)
        current = new_distance
        if numb == 0:
           min_distance = new_distance
        else:
           new_distance = min_distance + 0.2 * numb
        print new_distance, current
        output.write(' %2.1f'%(new_distance))
        dest_dir = workdir+'/'+str(new_distance)
        os.makedirs(dest_dir)
        numb += 1
        os.system('cp md.inp run-cp2k.sub '+dest_dir)
        os.system('mv geometry_1.xyz '+dest_dir+'/geometry.xyz')
if __name__ == '__main__':
    main()
    
