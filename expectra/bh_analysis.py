#!/usr/bin/env python
"""
This code is used to analyze the distribution of energy and s (chi_deviation).
How to use:
    gr_ES [pot_file] [traj_file] [de]
"""
import sys
import os
#from expectra.io import read_dots
from expectra.io import read_atoms
import numpy as np

def read_dots(filename):
    f = open(filename)
    energy = []
    accept=[]
    for line in f:
        if line.startswith('#'):
           continue
        fields = [ field for field in line.split()]
        energy.append(float(fields[6]))
        #s.append(float(fields[7]))
        accept.append(float(fields[2]))
        #dot = [float(fields[6]), float(fields[7])]
    f.close()
    return energy, accept

def log_atoms(energy, traj):
    log_trajectory = open('global_minimum.xyz', 'w')
    for i in range(len(energy)):
        log_trajectory.write("%d\n" % (len(traj[i])))
        log_trajectory.write("potential: %15.6f \n"
                                  %(energy[i]))
        for atom in traj[i]:
           log_trajectory.write("%s  %15.6f  %15.6f  %15.6f\n" % (atom.symbol,
                                    atom.x, atom.y, atom.z))
    log_trajectory.close()

def main():
    arg = sys.argv
    pot_file = arg[1]
    traj_file = arg[2]
    delt_e = float(arg[3])

    energy, accept = read_dots(pot_file)
    traj = read_atoms(traj_file)
    min_index = np.argmin(np.array(energy))
    min_energy = energy[min_index]

    #analyze
    numb_accept = 0.0
    numb_min = 0.0
    global_e = []
    global_traj =[]
    print len(energy), len(traj)
    new_min = False
    for i in xrange (len(energy)):
        numb_accept = numb_accept + accept[i]
        if energy[i] >= min_energy and energy[i] < min_energy + delt_e:
           numb_min = numb_min+1
           global_e.append(energy[i])
           global_traj.append(traj[i])
        if not new_min:
           if energy[i] >= min_energy and energy[i] < min_energy + 0.005:
              min_index = i
              new_min = True
        
    print "acceptace: ", numb_accept / float(len(accept)), min_index, min_energy, numb_min/float(len(energy))
    log_atoms(global_e, global_traj)

if __name__ == '__main__':
    main()

        
