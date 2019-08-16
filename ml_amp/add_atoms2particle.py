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
    atoms = read('POSCAR')
#    atoms.center()
    layers = find_layers(atoms.copy())
    inner_site = {}
    n_site = 0
    for key in layers.keys():
        new_sites, n_site = find_positions(layers[key][-1], n_site)
        if key==0:
           outer_site = new_sites
           continue
        inner_site.update(new_sites)
    n_Pd = 147
    n_H = int(n_Pd * 0.4)
    n_inners = int(n_H * 0.5)
    #print inner_site.keys()
    #inners = random.sample(random.shuffle(inner_site.keys()), n_inners)  
    #outers = random.sample(random.shuffle(outer_site.keys()), n_H - n_inners)
    inner_keys = inner_site.keys()
    outer_keys = outer_site.keys()
    #np.random.shuffle(inner_keys)
    #np.random.shuffle(outer_keys)
    #print inner_keys
    inners = random.sample(inner_keys, n_inners)  
    outers = random.sample(outer_keys, n_H - n_inners)
    print outers
    print inners
    #print n_inners
    for key in inners:
        if inner_site[key][0]=='bridge':
           atoms.append(Atom('H',inner_site[key][1]+inner_site[key][2] * 0.4))
        if inner_site[key][0]=='hollow':
           if inner_site[key][3]> 0:
               atoms.append(Atom('H',inner_site[key][1] + inner_site[key][2] * 0.6))
           else:
               atoms.append(Atom('H',inner_site[key][1] - inner_site[key][2] * 0.6))
    for key in outers:
        if outer_site[key][0]=='bridge':
           atoms.append(Atom('H',outer_site[key][1]+outer_site[key][2] * 0.4))
        if outer_site[key][0]=='hollow':
           if outer_site[key][3]>0:
               atoms.append(Atom('H',outer_site[key][1] + outer_site[key][2] * 0.5))
           else:
               atoms.append(Atom('H',outer_site[key][1] - outer_site[key][2] * 0.5))
    write('CONTCAR',atoms,format='vasp')
if __name__ == '__main__':
    main()

