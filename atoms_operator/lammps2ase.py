#!/usr/bin/env python
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.calculators.lammpsrun import Prism
from ase import Atoms
import numpy
import ase
import sys

def read_lammps_trj(filename=None, skip=0, every=1, maximum=1E36, specorder=None):
    """Method which reads a LAMMPS dump file."""
    if filename is None:
        print "No trajectory file is provided"

    if isinstance(specorder, str):
       specorder = specorder.split()
       
    atoms=[]
    f = open(filename, 'r')
    n_atoms = 0
    count = 0
    itrj = 0
    while True:
        line = f.readline()

        if not line or itrj > maximum:
            break

        if 'ITEM: TIMESTEP' in line:
            lo = [] ; hi = [] ; tilt = []
            id = [] ; type = []
            positions = [] ; velocities = [] ; forces = []
            # xph: add charges
            charges = []
            line = f.readline()
            #lei: itrj used for skipping
            itrj = int(line.split()[0])

        line = f.readline()
        if 'ITEM: NUMBER OF ATOMS' in line:
            line = f.readline()
            n_atoms = int(line.split()[0])
        
        #lei: skip geometries
        if itrj < skip:
           for i in range(n_atoms + 5):
               line = f.readline()
        else:
           line = f.readline()
           if 'ITEM: BOX BOUNDS' in line:
               # save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
               tilt_items = line.split()[3:]
               for i in range(3):
                   line = f.readline()
                   fields = line.split()
                   lo.append(float(fields[0]))
                   hi.append(float(fields[1]))
                   if (len(fields) >= 3):
                       tilt.append(float(fields[2]))
           
           line = f.readline()
           if 'ITEM: ATOMS' in line:
               # (reliably) identify values by labels behind "ITEM: ATOMS" - requires >=lammps-7Jul09
               # create corresponding index dictionary before iterating over atoms to (hopefully) speed up lookups...
               atom_attributes = {}
               for (i, x) in enumerate(line.split()[2:]):
                   atom_attributes[x] = i
               for n in range(n_atoms):
                   line = f.readline()
                   fields = line.split()
                   id.append( int(fields[atom_attributes['id']]) )
                   type.append( specorder[int(fields[atom_attributes['type']])-1] )
                   positions.append( [ float(fields[atom_attributes[x]]) for x in ['x', 'y', 'z'] ] )
#                   velocities.append( [ float(fields[atom_attributes[x]]) for x in ['vx', 'vy', 'vz'] ] )
#                   forces.append( [ float(fields[atom_attributes[x]]) for x in ['fx', 'fy', 'fz'] ] )
#                   if hasattr('charges'):
#                        charges.append(  float(fields[atom_attributes['q']]) )

               xhilo = (hi[0] - lo[0])
               yhilo = (hi[1] - lo[1])
               zhilo = (hi[2] - lo[2])
           
               cell = [[xhilo,0,0],[0,yhilo,0],[0,0,zhilo]]
           
               sort_atoms = sorted(zip(id, type, positions))
               
               cell_atoms = numpy.array(cell)
               type_atoms = [ types for (ids, types, position) in sort_atoms]
               positions = [ position for (ids, types, position) in sort_atoms]
               #forces = [ force for (ids, types, position) in sort_atoms]

               positions_atoms = numpy.array(positions)
               #forces_atoms = numpy.array(forces)
           
#               positions_atoms = np.array( [np.dot(np.array(r), rotation_lammps2ase) for r in positions] )
#               velocities_atoms = np.array( [np.dot(np.array(v), rotation_lammps2ase) for v in velocities] )
#               forces_atoms = np.array( [np.dot(np.array(f), rotation_lammps2ase) for f in forces] )

               count += 1
               if count % every == 0: 
                  atoms.append(Atoms(type_atoms, positions=positions_atoms, cell=cell_atoms, pbc=True))
    f.close()
    return atoms


def main():
   arg = sys.argv
   inputfile = arg[1]
   spec_one = arg[2]
   print [spec_one]
   spec_two = arg[3]
   outputfile = arg[4]
   p1 = read_lammps_trj(filename=inputfile,
                        skip = int(arg[5]),
                        every = 1,
                        maximum = int(arg[6]),
                        specorder = [spec_one, spec_two])
   cm = p1[0].get_center_of_mass()
   print cm
#   log_structures = Trajectory(filename+'.traj',
#                            'w', atoms)
   for p in p1:
      p.wrap(center=(0.5,0.5,0.5))
      #p.center(about=cm)  
   write(filename=outputfile,images=p1,format=outputfile.split('.')[1])
#   write(filename='au.xyz',images=p1,format='xyz')
if __name__ == '__main__':
    main()
