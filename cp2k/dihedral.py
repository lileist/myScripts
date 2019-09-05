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
    da = float(paras['step_size'])
    angle_span = 720
    trajs = read(filename=arg[2], index=':', format='xyz')
    output = open('w_distribution.dat','w')
    evolve_t = open('evolve_t.dat','w')
    #initialize
    gr_1 = []
    gr_2 = []
    #gr_2d = [[0 for x in range(w)] for y in range (h)]
    image_numb = len(trajs)
    t_numb_0 = 0
    i=0
    dihedrals = {}
    grs = {}
    t_numb = {}
    for target in targets:
        dihedrals[target] = []
        t_numb[target] = 0
        grs[target] = []

    for traj in trajs:
        #f_p1_v = traj[f_p1].position
        #norm_vect = numpy.cross(traj[f_p2].position - f_p1_v, traj[f_p3].position - f_p1_v)
        for target in targets:
           dihedral = get_dihedral_angle(traj, f_p1, f_p2, f_p3, target)
           dihedrals[target].append(dihedral)
    #                      * numpy.sign(numpy.dot(norm_vect, traj[target_1].position - f_p1_v)))
#    for i in xrange (image_numb):
#        for k in xrange (h):
#            for j in xrange (w):
#                if dihedral_1[i] >= angle_min + float(j)*da and dihedral_1[i] < angle_min + float((j+1))*da \
#                   and dihedral_2[i] >= angle_min + float(k)*da and dihedral_2[i] < angle_min+float((k+1))*da:
#                   gr_2d[k][j] += 1
 #                  t_numb_0 += 1
    #normalize
    maxs = -400.0
    mins = 400.0
    for target in targets:
       temp = numpy.array(dihedrals[target])
       #dihedrals[target] = temp - numpy.mean(temp)
       dihedrals[target] = temp
       temp_max = numpy.amax(dihedrals[target]) 
       temp_min = numpy.amin(dihedrals[target])
       if temp_max > maxs:
          maxs = temp_max
       if temp_min < mins:
          mins = temp_min
    
    w = int((maxs-mins)/da)+1
    h = w
    entropies = {}
    deviations = ''
    head = ''
    for target in targets:
      entropies[target]=0
      for i in range(w):
          grs[target].append(0.0)
      deviations += " {:8.4f} ".format(numpy.std(dihedrals[target]))
      head += " {:>8s} ".format(trajs[0][target].symbol+str(target))

    for i in range(image_numb):
        output_dihedrals = " {:8.4f} ".format(i*0.0005)
        for target in targets:
           dihedral = dihedrals[target][i]
           for j in range(w):
               if dihedral >= mins + float(j)*da and dihedral < mins + float((j+1))*da:
                  grs[target][j] += 1
                  t_numb[target] += 1
           output_dihedrals += " {:8.4f} ".format(dihedral)
        evolve_t.write("%s\n"% (output_dihedrals))
        evolve_t.flush()

    for i in range(w):
       angle = mins + float(i)*da 
       output_grs = " {:15.6f} ".format(angle)
       for target in targets:
          grs[target][i] = float(grs[target][i])/float(t_numb[target])
          output_grs += " {:15.6f} ".format(grs[target][i])
          if grs[target][i] > 0:
            entropies[target] += grs[target][i] * (-math.log10(grs[target][i]))
          #elif grs[target][i] < 0:
          #  entropy += -grs[target][i] * (-math.log10(-grs[target][i]))
       output.write('%s\n' % (output_grs))
#       for j in range(h):
#          gr_2d[i][j] = float(gr_2d[i][j])/float(t_numb)
    max_s = 0
    for target in targets:
       head += " {:>8s} ".format(trajs[0][target].symbol+str(target)+'_s')
       deviations += " {:8.4f} ".format(entropies[target])
       if max_s < entropies[target]:
          max_s = entropies[target]
    head += " {:>8s} ".format('max_s')
    deviations += " {:8.4f} ".format(max_s)
    log_stat = open('entropy.dat','w')
    log_stat.write("%s\n"% (head))
    log_stat.write("%s\n"% (deviations))
#    a=numpy.array(gr_2d)
#    printMatrix(a, angle_min, angle_max, angle_min, angle_max)

if __name__ == '__main__':
    main()
    
