"""
This code is used to adde H2 molecule on a particle
"""
import ConfigParser
from ase.io import read, Trajectory
from ase.optimize.fire import FIRE
from ase.atom import Atom
import numpy as np
import ConfigParser

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
    h_index=[]
    for index in range(len(atoms[0])):
        if atoms[0][index].symbol=='H':
            h_index.append(index)
    h_index.sort(reverse=True)
    for p in atoms:
        for index in h_index:
           p.pop(i=index) 
        center = np.mean(p.get_positions(),axis=0)
        print center
        for atom in p:
           distance = np.linalg.norm(center - atom.position)
           if distance > 2.0:
              dot = atom.position + bondlength * (atom.position-center)/distance
              perp_vector = perpendicular_vector(atom.position-center)
              unit_perp_vector = perp_vector/np.linalg.norm(perp_vector)
              for h_distance in h_distances:
                 image = p.copy()
                 image.append(Atom('H',dot + h_distance*unit_perp_vector / 2.0))
                 image.append(Atom('H',dot - h_distance*unit_perp_vector / 2.0))
                 #image=image.wrap(center=(0.5, 0.5, 0.5))
                 traj.write(image)
              break

if __name__ == '__main__':
    main()

