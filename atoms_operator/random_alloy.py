import numpy as np
from ase.build import *
from ase.constraints import FixAtoms
from ase.io import write
from ase.build import sort

seed=3021323
host='Ag'
impurity='Pt'

latticeconstant=4.0853*0.5+3.80*0.5
concentration_impurity=0.5


model=fcc111(host, size=(3,3,4), a=latticeconstant, vacuum=6.0, orthogonal=False)
c = FixAtoms(mask=[x >2   for x in model.get_tags()])
model.set_constraint(c)

elements=model.get_chemical_symbols()
num_atom=model.get_number_of_atoms()

num_impurity=np.round(num_atom*concentration_impurity)


np.random.seed(seed)

i=0
while i < int(num_impurity):
    r=np.random.rand()
    n=int(np.round(r*num_atom))
    if elements[n]==host:
        elements[n]=impurity
        i=i+1


model.set_chemical_symbols(elements)
model=sort(model)
write('POSCAR',model,format='vasp',direct=True)
