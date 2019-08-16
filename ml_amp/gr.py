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
    while True:
        elements = []
        positions = []
        line = f.readline()
        if not line:
           break
        if cycle == -1:
           atom_numb = int(line.split()[0])
        line = f.readline()
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
    try:
       configs = read(arg[1], index=":")
    except:
       configs = read_atoms(arg[1])

   #initialized based on chemical symbols
    dist_gr = {} 
    angle_gr = {}
    numb = {}
    numb_angle = {}
    output = {}
    output_angle = {}
    symbol_index={}
#    chem_symbols = configs[0].get_chemical_symbols()
    chem_symbols = ['H','Pd']
    for s_a in chem_symbols:
       for s_b in chem_symbols:
           key = [s_a,s_b]
           key = ''.join(sorted(key))
           dist_gr[key]=[]
           for s_c in chem_symbols:
               sides = sorted([s_a, s_c])
               key=''.join([sides[0], s_b, sides[1]])
               if key in angle_gr:
                  continue
               angle_gr[key]=[]
#    for config in configs[0]:
#        if config.symbol not in symbol_index:
#           symbol_index[config.symbol] = []
#           continue
#        symbol_index[config.symbol].append(config.index)
            
    dist_min = 0.0
    gr_numb = 400
    de = 0.02
    angle_min = 0.0
    gr_angle_numb = 360
    da = 1.0
    dist_record = open('dist.dat','w')
    for key in dist_gr:
       output[key] = open(key+'.dat', 'w')
       numb[key]=0
       for i in range(gr_numb):
           dist_gr[key].append(0.0)
    for key in angle_gr:
       output_angle[key] = open('angle_'+key+'.dat','w')  
       numb_angle[key] = 0
       for i in range(gr_angle_numb):
           angle_gr[key].append(0.0)
    natom = len(configs[0])
    nimage = 0
    t_dist = 0
    n_bond = 0
    Rc = 6.0

    short_PdPd = 0
    short_HPd = 0
    for config in configs:
        nimage += 1
        print nimage
        PdPd_checked = False
        HPd_checked = False
        for i in range(natom-1):
            for j in range(i+1,natom):
                key = [config[i].symbol,config[j].symbol]
                key = ''.join(sorted(key))
                ij_dist = config.get_distance(i, j)
                t_dist += ij_dist
                n_bond += 1
                #Check number of structures with distance below certain value
                if ij_dist <= 2.4 and not PdPd_checked:
                   short_PdPd+=1
                   PdPd_checked = True
                if ij_dist <= 1.4 and not HPd_checked:
                   short_HPd+=1
                   HPd_checked = True

                for k in range (gr_numb):
                    if ij_dist >= dist_min + float(k)*de and ij_dist < dist_min + float((k+1))*de:
                       dist_gr[key][k] += 1
                       numb[key]=numb[key]+1
                       break
  
    for key in dist_gr:
       for i in range(gr_numb):
           dist_gr[key][i] = float(dist_gr[key][i])/numb[key]
           r = dist_min + float(i)*de - de*0.5
           output[key].write('%15.6f  %15.6f\n' % (r, dist_gr[key][i]))
    print short_PdPd, short_HPd

if __name__ == '__main__':
    main()
    
