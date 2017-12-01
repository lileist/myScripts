#!/usr/bin/env python
"""
Code used to analyze DOSCAR and output pdos for the selected atoms
pdos.py [atoms_selected]
Files required: POSCAR DOSCAR PROCAR
Two ways to define the selected atoms:
    1. provide the index list of the atoms you're interested:
       pdos.py 121 122
    2. provide the atom type and a box which is defined by minimum and maximum xyz.
       pdos.py Ti 
       The code will ask you to input minimum and maximum xyz values to define the box
"""
from ase import Atoms
from ase.io import read
import sys
import os
import numpy as np
import copy

def pick_atoms(atoms, atom_type, 
               x_max = 100.0, 
               x_min = -100.0, 
               y_max = 100.0, 
               y_min = -100.0, 
               z_max = 100.0, 
               z_min = -100.0):
    selected_atoms = []
    for atom in atoms:
       if atom.symbol in atom_type:
          if atom.x < x_max and atom.x > x_min:
             if atom.y < y_max and atom.y > y_min:
                if atom.z < z_max and atom.z > z_min:
                   #print atom.index
                   selected_atoms.append(atom.index)
    return selected_atoms

    """
    f= open('POSCAR', 'r')
    lattice = []
    f.readline()
    f.readline()
    for i in range (3):
        line = f.readline()
        fields = [float(field) for field in line.split()]
        lattice.append(fields[i])
    line = f.readline()
    elements = [field for field in line.split()]
    line = f.readline()
    elem_number = [float(field) for field in line.split()]
    t_number = sum (elem_number)
    f.readline()
    if f.readline() == 'Direct':
       x_max = x_max / lattice[0]
       x_min = x_min / lattice[0]
       y_max = y_max / lattice[1]
       y_min = y_min / lattice[1]
       z_max = z_max / lattice[2]
       z_min = z_min / lattice[2]
    for i in range (t_number):
       line = f.readline()
       fields = [float(field[i]) for i in len(line.split()) if i < 3]
       if fields[0] 
    """        
    
def read_doscar(filename, selected_atoms=[]):
    f = open('DOSCAR', 'r')
    fprocar = open('PROCAR', 'r')
    #pdos_file = open(filename, 'w')
    tdos_file = open(filename, 'w')
    atom_numb = None
    for i in range (5): 
        f.readline()
        fprocar.readline()
    #read procar
    fprocar.readline()
    fprocar.readline()
    line = fprocar.readline()
    orbitals = [ field for field in line.split()]
    print orbitals
    pdos_up = {} 
    pdos_down = {} 
    n_l = len(orbitals) - 2
    cycle = -1
    initial_pdos = []
    while True:
        line = f.readline()
        if not line:
           break
        fields = [ float(field) for field in line.split()]
        nedos = int(fields[2])
        fermi = fields[3]
        #read and write total dos
        if cycle == -1:
           energy = []
           print "nedos:", nedos
           print "fermi:", fermi
           for i in range (nedos):
               line = f.readline()
               fields = [ float(field) for field in line.split()]
               tdos_file.write("%8.3f %15.8f %15.8f\n" %
                               (fields[0]-fermi, fields[1], fields[2]))
               energy.append(fields[0]-fermi)
               initial_pdos.append(0.0)
           #initialize pdos
           for i in range (n_l):
               pdos_up[orbitals[i+1]] = [2*i+1, np.array(initial_pdos)]
               pdos_down[orbitals[i+1]] = [2*(i+1), np.array(initial_pdos)]
           print pdos_up.keys()
           print pdos_up['s'][0], pdos_up['py'][0], pdos_up['pz'][0], pdos_up['px'][0], pdos_up['dxy'][0], pdos_up['dyz'][0], pdos_up['dz2'][0], pdos_up['dxz'][0], pdos_up['dx2'][0]
           print pdos_down['s'][0], pdos_down['py'][0], pdos_down['pz'][0], pdos_down['px'][0], pdos_down['dxy'][0], pdos_down['dyz'][0], pdos_down['dz2'][0], pdos_down['dxz'][0], pdos_down['dx2'][0]
        #read and store pdos for selected atoms
        else:
           if cycle in selected_atoms:
              #print 'n_l', n_l
              for i in range (nedos):
                  line = f.readline()
                  fields = [ float(field) for field in line.split()]
                  for key in pdos_up:
                     pdos_up[key][1][i] += fields[pdos_up[key][0]]
                     pdos_down[key][1][i] += fields[pdos_down[key][0]]
           else:
              for i in range (nedos):
                  f.readline()
        cycle += 1
    return energy, fermi, pdos_up, pdos_down

def assign_min_max(inputs):
    fields = [ float(field) for field in inputs.split()]
    if len(fields) > 0:
       min = fields[0]
       max = fields[1]
    else:
       min = -1000.0
       max = 1000.0
    return min, max
           
def main():
    arg = sys.argv
    atom_type = []
    for i in range (len(arg)):
        if i > 0 :
           atom_type.append(arg[i]) 
    index = True
    try:
       #atoms selected based on provided index
       selected_atoms = [int(field) for field in atom_type]
    except:
       index = False
       pass
    #atoms selected based on defined box and atom_type
    vacuum_level = raw_input("please enter vaccum level energy\n")
    vacuum_level = float(vacuum_level)
    if not index:
       x = raw_input("please enter x min and max\n")
       y = raw_input("please enter y min and max\n")
       z = raw_input("please enter z min and max\n")
       x_min, x_max = assign_min_max(x)
       y_min, y_max = assign_min_max(y)
       z_min, z_max = assign_min_max(z)
       print x_min, x_max, y_min, y_max, z_min, z_max
       p1 = read(filename = 'POSCAR', index=0, format = 'vasp')
       selected_atoms = pick_atoms(p1, atom_type, x_max, x_min, y_max, y_min, z_max, z_min)
    print selected_atoms
    energy, fermi, pdos_up, pdos_down = read_doscar('t_dos.dat', selected_atoms)

    pdosfile = open( ''.join(atom_type) + '_pdos.dat', 'w')
    print pdosfile
    s_dos_up   = []
    s_dos_down = []
    p_dos_up   = []
    p_dos_down = []
    d_dos_up   = []
    d_dos_down = []
    pdosfile.write("#%-8s %15s %15s %15s %15s %15s %15s\n" %
                    ('Energy', 's_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down'))
    for i in range (len(energy)):
        s_dos_up.append( pdos_up['s'][1][i] )
        s_dos_down.append( pdos_down['s'][1][i] )
        p_dos_up.append( pdos_up['py'][1][i] + pdos_up['pz'][1][i] + pdos_up['px'][1][i] )
        p_dos_down.append( pdos_down['py'][1][i] + pdos_down['pz'][1][i] + pdos_down['px'][1][i] )
        d_dos_up.append( pdos_up['dxy'][1][i] + pdos_up['dyz'][1][i] +pdos_up['dz2'][1][i] 
                          + pdos_up['dxz'][1][i] + pdos_up['dx2'][1][i] )
        d_dos_down.append( pdos_down['dxy'][1][i] + pdos_down['dyz'][1][i] +pdos_down['dz2'][1][i] 
                          + pdos_down['dxz'][1][i] + pdos_down['dx2'][1][i])
        pdosfile.write("%8.3f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n" % 
                        (energy[i]+fermi-vacuum_level, s_dos_up[i], -s_dos_down[i], p_dos_up[i], 
                         -p_dos_down[i], d_dos_up[i], -d_dos_down[i], 
                         s_dos_up[i]+ p_dos_up[i]+d_dos_up[i],
                         -s_dos_down[i]-p_dos_down[i]-d_dos_down[i])) 
    s_dos = np.add(np.array(s_dos_up), np.array(s_dos_down))
    p_dos = np.add(np.array(p_dos_up), np.array(p_dos_down))
    d_dos = np.add(np.array(d_dos_up), np.array(d_dos_down))
    s_center = np.dot(s_dos, np.array(energy)) / sum(s_dos)
    p_center = np.dot(p_dos, np.array(energy)) / sum(p_dos)
    if sum(d_dos) == 0.0:
       d_center = 0.0
    else:
       d_center = np.dot(d_dos, np.array(energy)) / sum(d_dos)
    pdosfile.write("%s %8.3f %8.3f %8.3f %8.3f\n" % ("#", fermi-vacuum_level, s_center, p_center, d_center))

    pdosfile.close()
if __name__ == '__main__':
    main()
