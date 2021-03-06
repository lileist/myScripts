from ase.lattice.surface import *
from ase.constraints import FixAtoms
from ase.io import write

#a=3.71
#slab = fcc111('Cu', size=(4,4,4), a=a, vacuum=10.0, orthogonal=True)
#slab = fcc100('Cu', size=(3,3,3), a=a, vacuum=10.0)
#slab = fcc211('Cu', size=(3,4,3), a=a, vacuum=5.0, orthogonal=True)
#a= 4.046
#slab = fcc111('Al', size=(4,6,6), a=a, vacuum=10.0, orthogonal=True)
#slab = fcc100('Cu', size=(3,3,3), a=a, vacuum=10.0)
#slab = fcc211('Cu', size=(3,4,3), a=a, vacuum=5.0, orthogonal=True)
#a= 3.92
#slab = fcc111('Pt', size=(4,4,4), a=a, vacuum=7.5, orthogonal=True)
#slab = fcc211('Pt', size=(6,5,3), a=a, vacuum=5.0, orthogonal=True)
a=4.08
#slab = fcc111('Au', size=(3,4,3), a=a, vacuum=7.5, orthogonal=True)
#slab = fcc111('Au', size=(4,4,3), a=a, vacuum=7.5, orthogonal=True)
slab = fcc211('Au', size=(6,4,3), a=a, vacuum=5.0, orthogonal=True)
#Ni
#a= 3.524
#slab = fcc111('Ni', size=(4,6,6), a=a, vacuum=20.0, orthogonal=True)
c = FixAtoms(mask=[x >2   for x in slab.get_tags()])
slab.set_constraint(c)
write('POSCAR',slab,format='vasp',direct=True)
