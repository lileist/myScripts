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

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    configs = read_atoms(paras['configs_file'])
    geo_center_1 = [int(field) for field in paras['geo_center_1'].split()]
    geo_center_2 = [int(field) for field in paras['geo_center_2'].split()]
    image_numb = 0
    output = open(paras['output_file'], 'w')
    #group atom: H2O, NH3
    for atoms in configs:
        geo_1 = numpy.array([0.,0.,0.])
        geo_2 = numpy.array([0.,0.,0.])
        for atom in atoms:
            if atom.index+1 in geo_center_1:
               geo_1 += atom.position
               print 'geo_1:',atom.index
            if atom.index+1 in geo_center_2:
               geo_2 += atom.position
               print 'geo_2:',atom.index
        geo_1 = geo_1 / float(len(geo_center_1))
        geo_2 = geo_2 / float(len(geo_center_2))
        print geo_1
        print geo_2
        output.write('%d %9.6f\n'%(image_numb,numpy.linalg.norm(geo_1 - geo_2) ))
        image_numb +=1
if __name__ == '__main__':
    main()
    
