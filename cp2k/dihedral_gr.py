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
    #if angle > 300:
    #   return -360+angle
    #else:
    return angle

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    targets = [int(field) for field in paras['target_atoms'].split()]
    f_p1 = int(paras['facet_p1'])
    f_p2 = int(paras['facet_p2'])
    f_p3 = int(paras['facet_p3'])
    angle_min = float(paras['minimum_angle'])
    angle_max = float(paras['maximum_angle'])
    #da = float(paras['step_size'])
    da = 5.0
    angle_span = 720
    output_gr = open('gr_3d.dat','w')
    evolve_t = open('evolve_t.dat','r')
    #initialize
    image_numb = 0
    t_numb_0 = 0
    dihedrals = {}
    grs = {}
    t_numb = {}
    for target in targets:
        dihedrals[target] = []
        t_numb[target] = 0
        grs[target] = []
    
    for line in evolve_t.readlines():
        i = 1 
        fields = [float(field) for field in line.split()]
        for target in targets:
           dihedrals[target].append(fields[i])
           i+=1
        image_numb += 1

    mins = {}
    maxs = {}
    w = {}
    for target in targets:
       temp = numpy.array(dihedrals[target])
       #dihedrals[target] = temp - numpy.mean(temp)
       dihedrals[target] = temp
       maxs[target] = numpy.amax(dihedrals[target]) 
       mins[target] = numpy.amin(dihedrals[target])
       w[target] = int((maxs[target]-mins[target])/da) + 1
       output_gr.write('# %5d %15.6f  %15.6f\n'%(target, mins[target], maxs[target]))

    gr_2d =[ [[ 0.0 for l in xrange(w[targets[2]]) ] for k in xrange(w[targets[1]])  ] for j in xrange(w[targets[0]])]
              
    for i in xrange (image_numb):
        print i
        for j in xrange (w[targets[0]]):
            for k in xrange (w[targets[1]]):
                for l in xrange (w[targets[2]]):
                   if dihedrals[targets[0]][i] >= mins[targets[0]] + float(j)*da and dihedrals[targets[0]][i] < mins[targets[0]] + float((j+1))*da \
                      and dihedrals[targets[1]][i] >= mins[targets[1]] + float(k)*da and dihedrals[targets[1]][i] < mins[targets[1]] + float((k+1))*da \
                      and dihedrals[targets[2]][i] >= mins[targets[2]] + float(l)*da and dihedrals[targets[2]][i] < mins[targets[2]] + float((l+1))*da:
                      gr_2d[j][k][l] += 1.0
                      t_numb_0 += 1.0
    #normalize
    entropy = 0
    deviations = ''
    head = ''
    for j in xrange (w[targets[0]]):
        for k in xrange (w[targets[1]]):
            for l in xrange (w[targets[2]]):
              gr_2d[j][k][l] =  gr_2d[j][k][l]/t_numb_0
              output_gr.write('%15.6f\n'%(gr_2d[j][k][l]))
              if gr_2d[j][k][l] > 0:
                entropy += gr_2d[j][k][l] * (-math.log10(gr_2d[j][k][l]))
    log_stat = open('3d_entropy.dat','w')
    log_stat.write("%8.4f\n"% (entropy))
#    a=numpy.array(gr_2d)
#    printMatrix(a, angle_min, angle_max, angle_min, angle_max)

if __name__ == '__main__':
    main()
    
