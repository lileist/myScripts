from ase.lattice.surface import *
from ase.constraints import FixAtoms
from ase.io import write
from ase.build import bulk

#a=4.2067452381
#a=4.08
p1=bulk('Pd', 'fcc', a=3.89, cubic=True)
write('POSCAR',p1,format='vasp',direct=True)
