from ase.io import read,write

atoms = read(filename='opt_lammps',index='::',format='lammps-dump')
write(filename='lammps_test.traj',images=atoms,format='traj')
