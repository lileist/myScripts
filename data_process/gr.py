#!/usr/bin/env python
"""
This code is used to calculate strain of nanoparticles.
Method: stretch partilces along certain directions by dx and fix the boundary atoms, do optimization
strain_cal_cluster [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    configs_file = POSCAR
    format = vasp
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
from ase import Atoms, Atom
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
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

def get_geometry_center(atoms):
    n = float(len(atoms))
    x=0
    y=0
    z=0
    for atom in atoms:
        x += atom.x
        y += atom.y
        z += atom.z
    return x/n, y/n, z/n

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    atoms = read(paras['configs_file'], index=0, format=paras['format'])
    dist_gr = []
    au_gr = []
    pt_gr = []
    numb = 0
    au_numb = 0
    pt_numb = 0
    dist_min = float(paras['minimum_distance'])
    gr_numb = int(paras['gr_number'])
    de = float(paras['step_size'])
    image_numb = 0
    for i in range(gr_numb):
        dist_gr.append(0.0)
        au_gr.append(0.0)
        pt_gr.append(0.0)
    #for atoms in configs:
    x, y, z = get_geometry_center(atoms)
    ps_atom = Atom('C', (x, y, z))
    atoms.append(ps_atom)
    atom_numb = len(atoms)
    for i in range(atom_numb-1):
        ij_dist = atoms.get_distance(i, atom_numb-1)
        for k in range (gr_numb):
            if ij_dist >= dist_min + float(k)*de and ij_dist < dist_min + float((k+1))*de:
               dist_gr[k] += 1
               numb+=1
               if atoms[i].symbol == 'Au':
                  au_gr[k] += 1
                  au_numb += 1
               if atoms[i].symbol == 'Pt':
                  pt_gr[k] += 1
                  pt_numb += 1
               break
    print au_numb, pt_numb, numb
    output = open(paras['outputfile'], 'w')
    for i in range(gr_numb):
        dist_gr[i] = float(dist_gr[i])/numb
        au_gr[i] = float(au_gr[i])/numb
        pt_gr[i] = float(pt_gr[i])/numb
        r = dist_min + float(i)*de - de*0.5
        output.write('%15.6f  %15.6f %15.6f %15.6f\n' % (r, dist_gr[i],au_gr[i],pt_gr[i]))
if __name__ == '__main__':
    main()
    
