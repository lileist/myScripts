"""
"""
from ase.io import read, Trajectory, write
from ase.optimize.fire import FIRE
from ase.atom import Atom
import numpy as np
from ase.neighborlist import NeighborList
import random, copy
import sys
from operator import itemgetter

def expand_particle(atoms, expansions, target):
    com=atoms.get_center_of_mass()
    p=atoms.copy()
    images = []
    for expansion in expansions:
      atoms = p.copy()
      pos = atoms[target].position
      vect = pos - com
      vect_norm = np.linalg.norm(vect)
      vect=vect/vect_norm
      atoms[target].position = pos + expansion * vect_norm * vect
      images.append(atoms)
    return images

def main():
    args = sys.argv
    imgs = read(args[1], index=":")
    random.shuffle(imgs)
    imgs_selected_Pd = imgs[:200]
    random.shuffle(imgs)
    imgs_selected_H = imgs[:100]
    expansions = [0.1*(i+1) for i in range(8)]
    traj=Trajectory('expanded.traj','w')
    for img in imgs_selected_Pd:
       Pds = [atom.index for atom in img if atom.symbol == 'Pd']
       com = img.get_center_of_mass()
       d_atoms = []
       for Pd in Pds:
          d_atoms.append([np.linalg.norm(img[Pd].position - com), Pd])
       ordered_atoms = sorted(d_atoms,key=itemgetter(0), reverse=False)
       imags_expanded=expand_particle(img, expansions, random.choice(ordered_atoms[-5:])[1])
       for im in imags_expanded:
          traj.write(im)

    for img in imgs_selected_H:
       Hs = [atom.index for atom in img if atom.symbol == 'H']
       com = img.get_center_of_mass()
       d_Hs = []
       for H in Hs:
          d_Hs.append([np.linalg.norm(img[H].position - com), H])
       ordered_atoms = sorted(d_Hs,key=itemgetter(0), reverse=False)
       imags_expanded=expand_particle(img, expansions, ordered_atoms[-1][1])
       for im in imags_expanded:
          traj.write(im)
if __name__ == '__main__':
    main()

