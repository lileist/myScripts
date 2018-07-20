#!/usr/bin/env python
"""
This code is used to wrap all trajectories back to unit cell
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

def main():
    arg = sys.argv
    configs = read(arg[1],index=':',format='xyz')
    output = open('wraped.xyz', 'w')
    for atoms in configs:
        atoms.set_cell([[19.9683,0,0],[0,19.9683,0],[0,0,35.0]],scale_atoms=False)
        atoms.set_pbc((True, True, True))
        atoms.wrap()
        output.write('%d\n'%(len(atoms)))
        output.write('\n')
        for atom in atoms:
           output.write('%s %15.6f %15.6f  %15.6f\n' % (atom.symbol, atom.x, atom.y, atom.z))
if __name__ == '__main__':
    main()
    
