#!/usr/bin/env python
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.calculators.lammpsrun import Prism
import ase
import sys

def write_lammps_data(fileobj, atoms, specorder=[], force_skew=False, write_charge=False):
    """Method which writes atomic structure data to a LAMMPS data file."""
    if isinstance(fileobj, str):
#        f = paropen(fileobj, 'w')
        f = open(fileobj, 'w')
        close_file = True
    else:
        # Presume fileobj acts like a fileobj
        f = fileobj
        close_file = False

    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError('Can only write one configuration to a lammps data file!')
        atoms = atoms[0]

    f.write(f.name + ' (written by ASE) \n\n')

    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)
    f.write('%d \t atoms \n' % n_atoms)

    if specorder is None:
        # This way it is assured that LAMMPS atom types are always
        # assigned predictively according to the alphabetic order 
        species = sorted(list(set(symbols)))
    else:
        # To index elements in the LAMMPS data file (indices must
        # correspond to order in the potential file)
        species = specorder
    n_atom_types = len(species)
    f.write('%d  atom types\n' % n_atom_types)

    p = Prism(atoms.get_cell())
    print p
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism_str()

    f.write('0.0 %s  xlo xhi\n' % xhi)
    f.write('0.0 %s  ylo yhi\n' % yhi)
    f.write('0.0 %s  zlo zhi\n' % zhi)
    
    if force_skew or p.is_skewed():
        f.write('%s %s %s  xy xz yz\n' % (xy, xz, yz))
    f.write('\n\n')

    f.write('Atoms \n\n')
    # xph: add charge in data file
    if write_charge:
        for i, r in enumerate(map(p.pos_to_lammps_str,
                                  atoms.get_positions())):
            s = species.index(symbols[i]) + 1
            charge = atoms[i].charge
            f.write('%6d %3d %.4f %s %s %s\n' % ((i+1, s, charge)+tuple(r)))
    else:
        for i, r in enumerate(map(p.pos_to_lammps_str,
                                  atoms.get_positions())):
            s = species.index(symbols[i]) + 1
            f.write('%6d %3d %s %s %s\n' % ((i+1, s)+tuple(r)))
    
    if close_file:
        f.close()


def main():
   if ase.__version__ > '3.12.0':
      print "please use ase version <= 3.12.0"
      sys.exit()
   arg = sys.argv
   inputfile = arg[1]
   spec_one = arg[2]
   print [spec_one]
   spec_two = arg[3]
   outputfile = arg[4]
   p1 = read(filename=inputfile,index=0,format=inputfile.split('.')[1])
   p1.set_cell([[20.8021,0,0],[0,20.8021,0],[0,0,20.8021]],scale_atoms=False)
   p1.set_pbc((True, True, True))
   for atom in p1:
      if atom.symbol == 'H':
         atom.charge = 0.42379999
      if atom.symbol == 'O':
         atom.charge = -0.84759998
#   for p in p1:
#      p.wrap(center=(0.5,0.5,0.5))
      #p.center(about=cm)  
   write_lammps_data(fileobj=outputfile, atoms=p1, specorder=[spec_one, spec_two], force_skew=False, write_charge=True)
#   write_lammps_data(fileobj=outputfile, atoms=p1, specorder=None, force_skew=False, write_charge=False)
if __name__ == '__main__':
    main()
