#from ase.lattice.surface import *
from ase.constraints import FixAtoms
from ase.io import write
from ase.build import surface
from ase import Atoms

a=4.085  #Au
b=4.078  #Ag
a_r = 3.0 
b_r = 1.0

c = a * a_r / ( a_r + b_r) + b * b_r/ (a_r + b_r)
print c
#slab = fcc111('Ag', size=(1,1,4), a=c, vacuum=7.5, orthogonal=False)
#slab = fcc100('Cu', size=(3,3,3), a=c, vacuum=10.0)
#slab = fcc211('Cu', size=(3,4,3), a=a, vacuum=5.0, orthogonal=True)
#a= 3.92
#slab = fcc111('Pt', size=(4,4,4), a=a, vacuum=7.5, orthogonal=True)
#slab = fcc211('Pt', size=(6,5,3), a=a, vacuum=5.0, orthogonal=True)
#a=4.08
#slab = fcc111('Au', size=(3,4,3), a=a, vacuum=7.5, orthogonal=True)
#slab = fcc111('Au', size=(4,4,3), a=a, vacuum=7.5, orthogonal=True)
Ag3Au = Atoms('Au3Ag',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[c, c, c],
              pbc=True)
s3 = surface(Ag3Au, (1, 1, 1), 4)
s3.center(vacuum=7.5, axis=2)

Ag3Au.set_chemical_symbols('AuAgAu2')
#s4 = surface(Pt3Rh, (2, 1, 1), 9)
#s4.center(vacuum=10, axis=2)

c = FixAtoms(mask=[x >2   for x in s3.get_tags()])
s3.set_constraint(c)
write('POSCAR',s3,format='vasp',direct=True)
