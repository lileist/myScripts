"""
"""
import ConfigParser
from ase.io import read, Trajectory, write
from ase.optimize.fire import FIRE
from ase.atom import Atom
import numpy as np
import ConfigParser
from ase.neighborlist import NeighborList
import random, copy
import sys

def expand_particle(atoms, expansions):
    com=atoms.get_center_of_mass()
    p=atoms.copy()
    images = []
    for expansion in expansions:
      atoms = p.copy()
      for atom in atoms:
         pos = atom.position
         vect = pos - com
         vect_norm = np.linalg.norm(vect)
         vect=vect/vect_norm
         atom.position = pos + expansion * vect_norm * vect
      images.append(atoms)
    return images

def main():
    args = sys.argv
    imgs = read(args[1], index=":")
    random.shuffle(imgs)
    imgs_selected = imgs[:200]
    expansions = [0.05*(i+1) for i in range(8)]
    traj=Trajectory('expanded.traj','w')
    for img in imgs_selected:
       imags_expanded=expand_particle(img, expansions)
       for im in imags_expanded:
          traj.write(im)
if __name__ == '__main__':
    main()

