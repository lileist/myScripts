#!/usr/bin/env python
"""
This code is used to extract geometry and exafs for given energy and s (chi_deviation).
How to use:
   dot_selection  [energy] [s] [pl_steps] [bh_steps] [species_number]
"""
import sys
import os
from ase.neighborlist import NeighborList
from ase.io import read, write

def main():
    arg = sys.argv
    atoms = read(filename=arg[1], format='xyz')
    cutoff = []
    coord_n = []
    n_pt = 0
    n_au = 0
    for i in range(len(atoms)):
        cutoff.append(float(arg[2]))
    nl = NeighborList(cutoff, skin=0, bothways=True, self_interaction=False)
    nl.update(atoms)
    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)
        coord_n.append([i,len(indices)])
        #coord_n[i]=len(indices)
        #coord_n.append(len(indices))
        if len(indices)<=9:
           if atoms[i].symbol == 'Pt':
              n_pt+=1
           if atoms[i].symbol == 'Au':
              n_au+=1
           atoms[i].symbol = 'Cu'
        if len(indices)==10:
           atoms[i].symbol = 'Pd'
    write('new.xyz', images=atoms,format='xyz')
    print float(n_pt)/float(n_pt+n_au),float(n_au)/float(n_pt+n_au)
    #print sorted(coord_n, key=lambda coord: coord[1])
if __name__ == '__main__':
    main()

        
