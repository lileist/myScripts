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

def push_apart(positions, pushapart=1.5,keeplist=None):
    movea = np.zeros(np.shape(positions))
    alpha = 0.025
    for w in range(500):
        moved = 0
        movea = np.zeros(np.shape(positions))
        for i in range(len(positions)):
            for j in range(i+1,len(positions)):
                d = positions[i] - positions[j]
                magd = np.sqrt(np.vdot(d,d))
                if magd == 0:
                   print i, j
                if magd < pushapart:
                    moved += 1
                    vec = d/magd
                    movea[i] += alpha *vec
                    movea[j] -= alpha *vec
        positions += movea
        if moved == 0:
            break
    return positions

def pull_closer(positions, neighbor_info, pullcloser):
    movea = np.zeros(np.shape(positions))
    alpha = 0.025
    print pullcloser
    for w in range(500):
        moved = 0
        movea = np.zeros(np.shape(positions))

        for key in neighbor_info.keys():
            for i in range(len(neighbor_info[key])):
                d =  positions[key] - positions[neighbor_info[key][i]]
                magd = np.sqrt(np.vdot(d,d))
                if magd == 0:
                   print key, neighbor_info[key][i]
                if magd > pullcloser[key][i]:
                    moved += 1
                    vec = d/magd
                    movea[key] -= alpha *vec
                    movea[neighbor_info[key][i]] += alpha *vec
        positions += movea
        if moved == 0:
            break
    return positions

def main():
    imgs = read('train.traj',index='0::10', format='traj')
    natoms = len(imgs[0])
    cutoff = 3.3
    mindist = 2.7

    traj = Trajectory('traj.traj','w')
    
    for atoms in imgs:
       nl=NeighborList([cutoff/2]*len(atoms), self_interaction=False, bothways=False)
       nl.update(atoms)
       neighbor_info = {}
       n_pairs = 0
       for i in range(natoms):
          i_indices, i_offsets = nl.get_neighbors(i)
          neighbor_info[i] = i_indices
          n_pairs += len(i_indices)
          print len(i_indices)
       #randomvalues=np.random.random((n_pairs,)) + mindist
       randomvalues=np.random.normal(mindist, 0.5, (n_pairs,))
       
       pullcloser = {}
       pointer = 0
       for key in neighbor_info.keys():
          pullcloser[key]=randomvalues[pointer:len(neighbor_info[key])+pointer:1]
          pointer += len(neighbor_info[key])
       
       atoms.set_positions(pull_closer(atoms.get_positions(), neighbor_info, pullcloser))
       #write('CONTCAR',atoms,format='vasp')
       traj.write(atoms)
if __name__ == '__main__':
    main()

