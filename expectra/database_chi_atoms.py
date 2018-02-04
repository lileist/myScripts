#!/usr/bin/env python
import sys
sys.path.append("../../")
import os

#from ase.test import NotAvailable
from ase.io import read
from ase.units import kB, fs
from ase.optimize.lbfgs import LBFGS
from ase.calculators.emt import EMT
from ase.atoms import Atoms
import numpy
from expectra.cal_exafs import Expectra
#from basin_surface import BasinHopping
from expectra.basin import make_dir
#from tsase.calculators.lammps_ext import LAMMPS

def log_visited_configs(visited_configs):
    data_base=visited_configs
    log_database = open('new_au147_all.xyz', 'w')
    log_exafs = open('new_au147_exafs_all.dat','w')
    for state in data_base:
        config =data_base[state]
        numb_atoms=len(config[4])
        log_database.write("%d\n"%(numb_atoms))
        log_database.write("images: %s energy: %15.6f chi_differ: %15.6f\n"%(state, config[0], config[1]))
        for atom in config[4]:
            log_database.write("%s %15.6f %15.6f %15.6f\n"%(atom.symbol, atom.x, atom.y, atom.z))
        log_database.flush()
        log_exafs.write("images: %s\n"%(state))
        for j in xrange(len(config[5])):
            log_exafs.write("%12.7f  %12.7f\n"%(config[5][j], config[6][j]))
        log_exafs.flush()
    log_database.close()
    log_exafs.close()
#prefix 'old' used to destinguish configs in databased and new one
def read_atoms(config_filename, exafs_filename, prefix=None):
    f = open(config_filename, 'r')
    f_exafs = open(exafs_filename,'r')
    configs={}
    while True:
        elements = []
        positions = []
        line = f.readline()
        if not line:
           break
        atom_numb = int(line.split()[0])
        line = f.readline()
        fields=line.split()
        state = prefix+fields[1]
        e_pot = float(fields[3])
        chi = float(fields[5])
        for i in range (atom_numb):
            line = f.readline()
            fields = line.split()
            elements.append(fields[0])
            positions.append( [ float(fields[j+1]) for j in range(3) ] )
        elements = numpy.array(elements)
        positions = numpy.array(positions)
        p = Atoms(elements, positions=positions)
        p.set_cell([[80,0,0],[0,80,0],[0,0,80]],scale_atoms=False,fix=None)
        configs[state] = [e_pot, chi,[0.0],0, p]
    while True:
        k=[]
        chi=[]
        line = f_exafs.readline()
        if not line:
           break
        image = prefix+line.split()[1]
        for i in range(382):
            line = f_exafs.readline()
            fields=line.split()
            k.append(float(fields[0]))
            chi.append(float(fields[1]))
        configs[image].append(numpy.array(k))
        configs[image].append(numpy.array(chi))
    f.close()
    return configs.copy()

vc_configs_1 = read_atoms('all_structures.xyz', 'all_exafs.dat',prefix='old_')
vc_configs_2 = read_atoms('all_structures.xyz_1', 'all_exafs.dat_1',prefix='old_1_')
vc_configs_1.update(vc_configs_2)

log_visited_configs(vc_configs_1)
