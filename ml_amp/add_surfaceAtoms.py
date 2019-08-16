"""
This code is used to adde H2 molecule on a particle
"""
import ConfigParser
from ase.io import read, Trajectory
from ase.optimize.fire import FIRE
from ase.atom import Atom
import numpy as np
import ConfigParser, random
from ase.neighborlist import NeighborList

def perpendicular_vector(v):
    if v[0]==0 and v[1]==0:
        if iszero(v[2]==0):
            # v is Vector(0, 0, 0)
            raise ValueError('zero vector')

        # v is Vector(0, 0, v.z)
        return np.array([0, 1, 0])

    return np.array([-v[1], v[0], 0])


def main():
    #config = ConfigParser.SafeConfigParser()
    config = ConfigParser.ConfigParser()
    config.read('config.ini')
    #print config
    main_paras = dict(config.items('main'))
    #print main_paras
    atoms = read(main_paras['structurefile'],
                index=main_paras['structure_slice'], 
                format=main_paras['fileformat'])
    h_distances = np.array([0.6+float(i)/10.0 for i in range(10)])  
    bondlength = 1.6
    traj = Trajectory('traj.traj','w')
    h_indices=[]
    for index in range(len(atoms[0])):
        if atoms[0][index].symbol=='H':
            h_indices.append(index)
    h_indices.sort(reverse=True)
    for img in atoms:
        p = img.copy()
        for index in h_indices:
           p.pop(i=index) 
        nl=NeighborList([3.2/2]*len(p), self_interaction=False, bothways=True)
        nl.update(p)
        center = np.mean(p.get_positions(),axis=0)
        sites = {}
        n_site = 0
        for atom in p:
           distance = np.linalg.norm(center - atom.position)
           i_indices, i_offsets = nl.get_neighbors(atom.index)
           if distance > 5.0 and len(i_indices) < 9:
              dot = atom.position + bondlength * (atom.position-center)/distance
              perp_vector = perpendicular_vector(atom.position-center)
              unit_perp_vector = perp_vector/np.linalg.norm(perp_vector)
              sites[n_site] = [dot, unit_perp_vector]
              n_site += 1
        h_delete = []
        h_kept = []
        keys_toadd = random.sample(sites.keys(),3)
        for h_index in h_indices:
          delete = False
          for key in keys_toadd:
             if np.linalg.norm(img[h_index].position - sites[key][0]) < 3.0:
                h_delete.append(h_index)
                delete = True
                break
          if not delete:
             h_kept.append(h_index)
           #image=image.wrap(center=(0.5, 0.5, 0.5))
        #Delete same number of H atoms to ensure constant number of H in particle
        if len(h_delete) > 6:
           print "Too many h to delete", len(h_delete)
           continue
        else:
           h_delete.extend(random.sample(h_kept, 6-len(h_delete)))
           #print len(h_delete)
           del img[h_delete]

        n_H=0
        for a in img:
          if a.symbol=='H':
            n_H+=1
        print n_H
        for h_distance in h_distances:
           image = img.copy()
           for key in keys_toadd:
                image.append(Atom('H',sites[key][0] + h_distance * sites[key][1] / 2.0))
                image.append(Atom('H',sites[key][0] - h_distance * sites[key][1] / 2.0))
           traj.write(image)

if __name__ == '__main__':
    main()

