import numpy as np

def push_apart(positions, pushapart=1.5):
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
    #atoms
    return positions
from ase.io import read, write
atoms = read('CONTCAR', format='vasp')
atoms.set_positions(push_apart(atoms.get_positions(), pushapart=1.5))
write('CONTCAR_pushed',atoms, format='vasp')
