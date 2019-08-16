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

def perpendicular_vector(v):
    if v[0]==0 and v[1]==0:
        if iszero(v[2]==0):
            # v is Vector(0, 0, 0)
            raise ValueError('zero vector')

        # v is Vector(0, 0, v.z)
        return np.array([0, 1, 0])

    return np.array([-v[1], v[0], 0])

def find_layers(atoms):
    cutoff = 3.2

    nlayer = 0
    layers = {}
    while True:
       if len(atoms) == 0:
          break
       nl=NeighborList([cutoff/2]*len(atoms), self_interaction=False, bothways=True)
       nl.update(atoms)
       com=atoms.get_center_of_mass()
       core = []
       for i in range(len(atoms)):
          #first_layer (most outer layer)
          i_indices, i_offsets = nl.get_neighbors(i)
          if i==13 or i==15 or i==3:
             print i, len(i_indices)
          if len(i_indices) < 10:
             if nlayer not in layers.keys():
                layers[nlayer] = [i] 
             else:
                layers[nlayer].append(i)
          else:
            core.append(i)
       #layers[nlayer].sort(reverse=True) 
       shell = atoms.copy()
       del shell[core]
       del atoms[layers[nlayer]]
       write('shell_'+str(nlayer)+'.xyz',shell,format='xyz')
       write('core_'+str(nlayer)+'.xyz',atoms,format='xyz')
       layers[nlayer].append(shell.copy())
       print nlayer
       print "     ",layers[nlayer]
       nlayer += 1
    return layers

def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if (a_set & b_set): 
        return (a_set & b_set) 
    else: 
        return None 

def find_positions(atoms, n_site=0):
    
    nO=10
    nl=NeighborList([3.3/2]*len(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    com=atoms.get_center_of_mass()

    hollow = {}
    n_hollow = 0
    bridge ={}
    n_bridge = 0
    top = {}
    n_top = 0
    sites = {}
    for i in range(len(atoms)):
       i_indices, i_offsets = nl.get_neighbors(i)
       for j in i_indices:
           #if i < j:
           #   bpos = (atoms.positions[i]+atoms.positions[j]) / 2.0
           #   sites[n_site] = ['bridge',bpos, bpos-com]
           #   n_site += 1
           j_indices, j_offsets = nl.get_neighbors(j)
           c = common_member(i_indices, j_indices)
           if c is not None:
               for k in c:
                   k_indices, k_offsets = nl.get_neighbors(k)
                   if len(i_indices) < 10 and len(j_indices) < 10 and len(k_indices) < 10 and i<j and j<k:
                       #n=n+1
                       pos = (atoms.positions[i]+atoms.positions[j]+atoms.positions[k])/3.
                       v1=atoms.positions[j]-atoms.positions[i]
                       v2=atoms.positions[k]-atoms.positions[i]
                       vn=np.cross(v1,v2)
                       vn=vn/np.linalg.norm(vn)
                       vo=pos-com
                       d=np.dot(vo,vn)
                       sites[n_site]=['hollow',pos, vn, d]
                       n_site += 1
    temp = copy.deepcopy(sites)
    """
    for i in range(n_site-1):
       for j in range(i+1,n_site):
          dist =np.linalg.norm(sites[i][1] - sites[j][1])
          if dist<0.5:
             print 'duplicated sites:', i, j, dist
             del temp[j]
    """
    return temp, n_site

def main():
    args = sys.argv
    imgs = read(args[1], index="::40")
    #layers = find_layers(atoms.copy())
    traj = Trajectory('traj.traj','w')

    H_indices = random.sample([a.index for a in imgs[0] if a.symbol == 'H'],8)

    n_img = 0
    for atoms in imgs:
       nl=NeighborList([2.5/2]*len(atoms), self_interaction=False, bothways=True)
       nl.update(atoms)
       pair_selected = []
       for H_index in H_indices:
         nl_indices, nl_offsets = nl.get_neighbors(H_index)
         pair_selected.append([H_index, random.sample(nl_indices, 1)[0]])
       for HPd_dist in [1.0, 1.1, 1.2, 1.3, 1.4]:
          img = atoms.copy()
          for pair in pair_selected:
            H_index = pair[0]
            Pd_selected = pair[1]
            v = atoms[H_index].position - atoms[Pd_selected].position
            vn = v/np.linalg.norm(v)
            del img[H_index]
            img.append(Atom('H',atoms[Pd_selected].position + vn * HPd_dist))
          traj.write(img)
          print n_img
          n_img+=1
if __name__ == '__main__':
    main()

