#!/usr/bin/env python
"""
This code is used to calculate the dihedral angle based on cp2k trajectory file
An example of inputfile (any line starts with '#' will be ignored):
 target_atom_1 = 0
 target_atom_1 = 1
 facet_p1 = 2    #define facet
 facet_p2 = 3
 facet_p3 = 4
 minimum_angle = -180
 maximum_angle = 180
 step_size = 1
"""
import sys
import os
import errno
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
from ase.constraints import constrained_indices, FixAtoms
import numpy, math
import subprocess

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
           raise

def readinputs(filename):
    f=open(filename, 'r')
    parameters = {}
    lines=f.readlines()
    for line in lines:
      if line.startswith('#'):
         continue
      fields = line.split('=')
      parameters[fields[0].strip()]=fields[1].replace("\n","").strip()
    return parameters

def printMatrix(a, e_min, e_max, s_min, s_max):
   print "# Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]"
   print "# "+ ("%15.6f" % e_min) +("  %15.6f" % e_max) +("  %15.6f" % s_min) +("  %15.6f" % s_max) 
#+("  %15.6f" % e_mean)+("  %15.6f" % s_mean)
   rows = a.shape[0]
   cols = a.shape[1]
   for i in range(0,rows):
      for j in range(0,cols):
         print("%8.5f" %a[i,j]),
      print
   print  

def get_dihedral_angle(atoms, f_p1, f_p2, f_p3, target_1):
    angle = atoms.get_dihedral(f_p1, f_p2, f_p3, target_1,mic=False)
    if angle > 300:
       return -360+angle
    else:
       return angle

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    target_1 = int(paras['target_atom_1'])
    target_2 = int(paras['target_atom_2'])
    f_p1 = int(paras['facet_p1'])
    f_p2 = int(paras['facet_p2'])
    f_p3 = int(paras['facet_p3'])
    angle_min = float(paras['minimum_angle'])
    angle_max = float(paras['maximum_angle'])
    da = float(paras['step_size'])
    angle_span = 360
    trajs = read(filename=arg[2], index=':', format='xyz')
    output = open('w_distribution.dat','w')
    evolve_t = open('evolve_t.dat','w')
    #initialize
    dihedral_1 = []
    dihedral_2 = []
    gr_1 = []
    gr_2 = []
    w = int(angle_span/da)+1
    h = w
    gr_2d = [[0 for x in range(w)] for y in range (h)]
    image_numb = len(trajs)
    t_numb=0
    t_numb_1 =0
    for i in range(w):
        gr_1.append(0.0)
    i=0
    for traj in trajs:
        f_p1_v = traj[f_p1].position
        norm_vect = numpy.cross(traj[f_p2].position - f_p1_v, traj[f_p3].position - f_p1_v)
        dihedral_1.append(get_dihedral_angle(traj, f_p1, f_p2, f_p3, target_1))
    #                      * numpy.sign(numpy.dot(norm_vect, traj[target_1].position - f_p1_v)))
        dihedral_2.append(get_dihedral_angle(traj,f_p1, f_p2, f_p3, target_2))
    #                      * numpy.sign(numpy.dot(norm_vect, traj[target_2].position - f_p1_v)))
        evolve_t.write("%8.4f  %8.4f %8.4f\n"% (i*0.0005, dihedral_1[i], dihedral_2[i]))
        for j in range(w):
            if dihedral_1[i] >= angle_min + float(j)*da and dihedral_1[i] < angle_min + float((j+1))*da:
               gr_1[j] += 1
               t_numb_1 +=1
        i+=1
        evolve_t.flush()
    for i in xrange (image_numb):
        for k in xrange (h):
            for j in xrange (w):
                if dihedral_1[i] >= angle_min + float(j)*da and dihedral_1[i] < angle_min + float((j+1))*da \
                   and dihedral_2[i] >= angle_min + float(k)*da and dihedral_2[i] < angle_min+float((k+1))*da:
                   gr_2d[k][j] += 1
                   t_numb+=1
    #normalize
    entropy = 0
    for i in range(w):
       gr_1[i] = float(gr_1[i])/float(t_numb_1)
       if gr_1[i] > 0:
         entropy += gr_1[i] * (-math.log10(gr_1[i]))
       angle = angle_min + float(i)*da 
       output.write('%15.6f  %15.6f\n' % (angle, gr_1[i]))
       for j in range(h):
          gr_2d[i][j] = float(gr_2d[i][j])/float(t_numb)
    deviation = numpy.std(dihedral_1)
    log_stat = open('entropy.dat','w')
    log_stat.write("%8.4f  %8.4f\n"% (entropy, deviation))
    a=numpy.array(gr_2d)
    printMatrix(a, angle_min, angle_max, angle_min, angle_max)

if __name__ == '__main__':
    main()
    
