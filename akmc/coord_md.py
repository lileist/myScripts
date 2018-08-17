"""
This code used to calculate coordination number of Au and number of Au atoms on the surface
mpirun -n 24 coord_md.py [trajectoryFilename] [skip] [every]
"""
#!/usr/bin/env python

import sys
import os
import re
import glob
import ast
#import pandas as pd
#import matplotlib.pyplot as plt
#from pele.storage import Database
import numpy
import mpi4py.MPI
from ase.neighborlist import neighbor_list as nl
from ase.io import read
from ase import Atoms
from operator import itemgetter
#from ase.calculators.lammpsrun import Prism
from lammpsrun import Prism
from tsase.calculators.lmplib import LAMMPSlib
from tsase.optimize.sdlbfgs import SDLBFGS

#write_lammps_data(filename, atoms, atom_types, comment=None, cutoff=None,
#                      molecule_ids=None, charges=None, units='metal')

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
            if not charge:
               charge = 0.0
            f.write('%6d %3d %.4f %s %s %s\n' % ((i+1, s, charge)+tuple(r)))
    else:
        for i, r in enumerate(map(p.pos_to_lammps_str,
                                  atoms.get_positions())):
            s = species.index(symbols[i]) + 1
            f.write('%6d %3d %s %s %s\n' % ((i+1, s)+tuple(r)))
    
    if close_file:
        f.close()

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
                   velocities.append( [ float(fields[atom_attributes[x]]) for x in ['vx', 'vy', 'vz'] ] )
                   forces.append( [ float(fields[atom_attributes[x]]) for x in ['fx', 'fy', 'fz'] ] )
#                   if hasattr('charges'):
#                        charges.append(  float(fields[atom_attributes['q']]) )

               xhilo = (hi[0] - lo[0])
               yhilo = (hi[1] - lo[1])
               zhilo = (hi[2] - lo[2])
           
               cell = [[xhilo,0,0],[0,yhilo,0],[0,0,zhilo]]
           
               sort_atoms = sorted(zip(id, type, positions, velocities, forces))
               
               cell_atoms = numpy.array(cell)
               type_atoms = [ types for (ids, types, position, velocity, force) in sort_atoms]
               positions = [ position for (ids, types, position, velocity, force) in sort_atoms]
               forces = [ force for (ids, types, position, velocity, force) in sort_atoms]

               positions_atoms = numpy.array(positions)
               forces_atoms = numpy.array(forces)
           
#               positions_atoms = np.array( [np.dot(np.array(r), rotation_lammps2ase) for r in positions] )
#               velocities_atoms = np.array( [np.dot(np.array(v), rotation_lammps2ase) for v in velocities] )
#               forces_atoms = np.array( [np.dot(np.array(f), rotation_lammps2ase) for f in forces] )

               count += 1
               if count % every == 0: 
                  atoms.append([itrj,Atoms(type_atoms, positions=positions_atoms, cell=cell_atoms, pbc=True)])
    f.close()
    return atoms

def geo_opt(atoms=None):
    data_file = 'data.lammps'
    write_lammps_data(fileobj=data_file,atoms=atoms, specorder=['Pd','Au'], force_skew=False, write_charge=True)
    os.system('lmp_serial <in.opt 1>out 2>err')
    atoms = read_lammps_trj(filename='trj_lammps', skip=0, every=1, specorder=['Pd','Au'])
#    cmds = ["pair_style eam/alloy",
#            "pair_coeff * * PdAu.set Pd Au"]
#    lammps = LAMMPSlib(lmpcmds = cmds, 
#                       atoms=atoms,
#                       lammps_header=['units metal',
#                                      'atom_style charge',
#                                      'atom_modify map array sort 0 0.0'])
#    atoms.set_calculator(lammps)
    #opt= SDLBFGS(atoms=atoms, maxmove = 1.0, dt = 0.2, dtmax = 1.0, logfile=None, trajectory=None)
#    opt= SDLBFGS(atoms=atoms, maxstep=0.2)
#    opt.run(fmax=0.02, steps=1000)
    return atoms[-1][1]

def get_coord(atoms):
    surf_Au = 0
    i = nl('i', atoms,
                     {('Au','Au'):3.4,
                     ('Au','Pd'):3.4,
                     ('Pd','Pd'):3.4
                     })
    coord = numpy.bincount(i)
    index_Au = [ atom.index for atom in atoms if atom.symbol=='Au']
    Au_coord=0
    for i in range(len(coord)):
       if i in index_Au:
          Au_coord += coord[i]
          if coord[i] < 10:
            surf_Au += 1
    return float(Au_coord)/float(len(index_Au)),surf_Au 

cwd = os.getcwd()

def bunch_cal_coord(configs, rank):
    results =[]
    if not os.path.exists(str(rank)):
       os.makedirs(str(rank))
    #have in.opt PdAu.set in directory
    os.system('cp in.opt PdAu.set '+str(rank))
    os.chdir(cwd+'/'+str(rank))
    for config in configs:
       print "Optimize Configs", config[0]
       Au_coord, surf_Au = get_coord(geo_opt(config[1]))
       results.append([config[0],Au_coord, surf_Au])
       print "   ",config[0],Au_coord, surf_Au
    os.chdir(cwd)
    return numpy.array(results)

comm = mpi4py.MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
root = 0
sendbuf = []

if rank == 0:
   args = sys.argv
   configs = read_lammps_trj(filename=args[1],
                        skip = int(args[2]),
                        every = int(args[3]),
                        maximum = int(args[4]),
                        specorder = ['Pd', 'Au'])
   #split data into chunks based on size of rank
   chunks = [[] for _ in range(size)]
   for i, chunk in enumerate(configs):
     chunks[i % size].append(chunk)
   sendbuf=chunks
#else:
#   sendbuf=None
#   v=None
v=comm.scatter(sendbuf, root)
v = bunch_cal_coord(v, comm.rank)
recvbuf=comm.gather(v,root)

if comm.rank==0:
   f = open("coord_time.dat",'w')
#   print len(recvbuf)
#   print recvbuf
#   recvbuf = sorted(recvbuf, key=itemgetter(0))
#   recvbuf = sorted(recvbuf, key=lambda x:x[0])
   for i in range(len(recvbuf)):
      for item in recvbuf[i]:     
         f.write("%20d  %6.4f  %6d\n"%(item[0], item[1], item[2]))
