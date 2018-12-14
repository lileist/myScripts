"""
This code used to calculate coordination number of Au and number of Au atoms on the surface
"""
#!/usr/bin/env python

import sys
import os
import re
import glob
import ast
import pandas as pd
import matplotlib.pyplot as plt
from pele.storage import Database
import numpy
import mpi4py.MPI
from ase.neighborlist import neighbor_list as nl
from ase.io import read
from ase import Atoms
from operator import itemgetter

md_time=[]
def read_lammps_trj(filename=None, skip=0, every=1, specorder=None):
    """Method which reads a LAMMPS dump file."""
    if filename is None:
        print "No trajectory file is provided"

    if isinstance(specorder, str):
       specorder = specorder.split()
       
    atoms=[]
    f = open(filename, 'r')
    n_atoms = 0
    count = 0
    while True:
        line = f.readline()

        if not line:
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
           md_time.append(itrj)
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

def get_coord(atoms):
    surf_Au = 0
    i = nl('i', atoms,
                     {('Au','Au'):3.3,
                     ('Au','Pd'):3.3,
                     ('Pd','Pd'):3.3
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

def bunch_cal_coord(configs):
    results =[]
    for config in configs:
       Au_coord, surf_Au = get_coord(config[1])
       results.append([config[0],Au_coord, surf_Au])
    return numpy.array(results)

args = sys.argv

configs = read_lammps_trj(filename=args[1],
                     skip = int(args[2]),
                     every = int(args[3]),
                     specorder = ['Pd', 'Au'])

comm = mpi4py.MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
root = 0
sendbuf = []

if rank == 0:
   #split data into chunks based on size of rank
   chunks = [[] for _ in range(size)]
   for i, chunk in enumerate(configs):
     chunks[i % size].append(chunk)
#   print("chunks:",chunks)
   sendbuf=chunks
#else:
#   sendbuf=None
#   v=None
v=comm.scatter(sendbuf, root)
v = bunch_cal_coord(v)
recvbuf=comm.gather(v,root)

if comm.rank==0:
   f = open("coord_time.dat",'w')
   print recvbuf
#   recvbuf = sorted(recvbuf, key=itemgetter(0))
#   recvbuf = sorted(recvbuf, key=lambda x:x[0])
   for i in range(len(recvbuf)):
      for item in recvbuf[i]:     
         f.write("%20d  %6.4f  %6d\n"%(item[0], item[1], item[2]))
