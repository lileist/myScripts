#!/usr/bin/env python
"""
This code is used to extract geometry, chi and energy information from expectra outputs
After exepctra runs, geometries stored in 'configs', chi in 'exafs', 'state, energy, S' in 'visited_configs
"""
import sys
import os
import copy
from expectra.io import read_chi 
import numpy 
from ase.atoms import Atoms
from expectra.atoms_operator import match, single_atom
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
def read_chi(filename):
    f = open(filename)
    ks = []
    chis = []
    for line in f:
        if "state" in line:
            continue
        fields = [ float(field) for field in line.split() ]
        k = fields[0]
        chi = fields[1]
        ks.append(k)
        chis.append(chi)
    f.close()
    ks = numpy.array(ks)
    chis = numpy.array(chis)
    return ks, chis

def read_atoms(filename):
    f = open(filename, 'r')
    elements = []
    positions = []
    while True:
        line = f.readline()
        if not line:
           break
        atom_numb = int(line.split()[0])
        line = f.readline()
        fields=line.split()
        e_pot = float(fields[3])
        for i in range (atom_numb):
            line = f.readline()
            fields = line.split()
            elements.append(fields[0])
            positions.append( [ float(fields[j+1]) for j in range(3) ] )
        elements = numpy.array(elements)
        positions = numpy.array(positions)
    p = Atoms(elements, positions=positions)
    p.set_cell([[80,0,0],[0,80,0],[0,0,80]],scale_atoms=False,fix=None)
    f.close()
    return p

def read_vc(directory):
    f = open(directory+'visited_configs.dat', 'r')
    print directory+'visited_configs.dat'
    configs=[]
    visited_configs={}
    while True:
        line = f.readline()
        fields = line.split()
        if 'pareto_step' in line:
           pl_step = int(fields[1])
           configs.append(visited_configs)
           visited_configs={}
           continue
        if not line:
           break
        #visited_configs[ fields[0] ]= [float(fields[1]), float(fields[2])]
        visited_configs[ fields[0] ] =[ float(fields[1]), float(fields[2])]
    latest_configs = []
    for state in configs[-1]:
        latest_configs.append([configs[-1][state][0], configs[-1][state][1],
                              read_atoms(directory+'configs/'+state), read_chi(directory+'exafs/'+state+'_Au')])
    f.close()
    return copy.deepcopy(latest_configs)

cwd = os.getcwd()
initial_direct = cwd+'/0/'
print initial_direct
vc_configs = read_vc(initial_direct)
print len(vc_configs)
data_direct = ['0.2', '0.4','0.6','0.8','1.0']
database = copy.deepcopy(vc_configs)

#state: [energy, chi, atoms, ks_array, chis_array]
for direct in data_direct:
    new_configs = read_vc(cwd+'/'+direct+'/')
    print len(new_configs)
    print cwd+'/'+direct+'/'
    for config in new_configs:
        matched = False
        for config_ref in vc_configs:
            if config[0] - config_ref[0] < 0.0001:
               if match(config[2], config_ref[2], 0.2, 3.0, True):
                  matched = True
                  break
        if matched:
           continue
        else:
           print "new_geometry found"
           database.append(config)

log_database = open(cwd+'/'+'au55_database_1.xyz', 'w')
log_exafs = open(cwd+'/'+'au55_exafs_1.dat','w')
numb_atoms=len(database[0][2])
for config in database:
    log_database.write("%d\n"%(numb_atoms))
    log_database.write("images: %d energy: %15.6f chi_differ: %15.6f\n"%(config.index(), config[0], config[1]))
    for atom in config[2]:
        log_database.write("%s %15.6f %15.6f %15.6f\n"%(atom.symbol(), atom.x, atom.y, atom.z))
    log_database.flush()
    log_exafs.write("images: %d\n"%(config.index()))
    for i in range(len(config[3])):
        log_exafs.write("%12.7f  %12.7f\n"%(config[3][i], config[4][i]))
    log_exafs.flush()
log_database.close()
        
           
           
               
    


