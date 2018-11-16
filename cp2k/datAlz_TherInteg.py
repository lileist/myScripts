#!/usr/bin/env python
"""
This code is used to calculate strain of nanoparticles.
Method: stretch partilces along certain directions by dx and fix the boundary atoms, do optimization
strain_cal_cluster [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    cp2k_inp = suppl.inp
    job_submit_script = qsub.hf
    job_submit_cmd    = sbatch
"""
import sys
import os
import errno
import ase
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp
from ase.utils.geometry import sort
#from ase.constraints import constrained_indices, FixAtoms
from numpy import linalg as la
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

def simpson(f, a, b, n):
    if n%2:
       raise ValueError("n must be even (received n=%d)" % n)
    h = (b - a) / n
    s= f[str(a)] +f[str(b)]
    
    for i in range(1, n, 2): 
        s += 4 * f[str(a + i * h)]
    for i in range(2, n-1, 2):
        s += 2 * f[str(a + i * h)]
    return s * h / 3
    
def find_distance(filename, index_1, index_2,max_numb):
    f = open(filename,'r')
    line_numb = 0
    distance = []
    while True:
       line = f.readline()
       if not line:
          break
       line_numb+=1
       if line_numb == 1:
          atom_numb = int(line.split()[0])
       line = f.readline()
       positions = []
       for i in range(atom_numb):
           line = f.readline()
           if i == index_1 or i == index_2:
              fields = line.split()
              positions.append( numpy.array([ float(fields[j+1]) for j in range(3) ] ))
       #distance.append(la.norm(positions[0]-positions[1]))
       distance.append(math.sqrt((positions[0][0]-positions[1][0])**2+(positions[0][1]-positions[1][1])**2+(positions[0][2]-positions[1][2])**2))
       if line_numb > max_numb:
          break
    print 'images',line_numb
    f.close()
    return numpy.array(distance)

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    distances=paras['distance'].split()
    output = open('freeEnergy.dat','w')
    dx = float(paras['dx'])
    avg_force = []
    free_energy = []
    md_steps = []
    atom_distance = {}
    std_data = open('avg_std.dat', 'w')
    force_file = None
    for distance in distances:
        if force_file is None:
           cp2k_inp = open(distance+'/'+paras['cp2k_inp'], 'r')
           lines = cp2k_inp.readlines()
           for line in lines:
               if 'PROJECT_NAME' in line:
                  force_file = line.split()[1]+'-1.LagrangeMultLog'
                  break
        input_force = open(distance+'/'+force_file, 'r')
        total_force = 0
        shake_numb = 0
        numb_force = 0
        numb_line = 0
        lines = input_force.readlines()
        #find distance between two atoms
        if len(arg) > 3:
           atom_distance[distance] = find_distance(distance+'/'+arg[3], int(paras['index_1']), int(paras['index_2']), int(paras['max_step']))
           if int(paras['max_step'])>= len(atom_distance[distance]):
              std_data.write("%s %15.6f %15.6f\n"%(distance, numpy.average(atom_distance[distance][int(paras['stablize_step']):]),numpy.std(atom_distance[distance][int(paras['stablize_step']):])))
           else:
              std_data.write("%s %15.6f %15.6f\n"%(distance, numpy.average(atom_distance[distance][int(paras['stablize_step']):int(paras['max_step'])]),numpy.std(atom_distance[distance][int(paras['stablize_step']):int(paras['max_step'])])))
        for line in lines:
#            numb_line += 1
            if 'Shake  Lagrangian' in line:
               shake_numb += 1
               if shake_numb < int(paras['stablize_step']) or shake_numb > int(paras['max_step']):
                  continue
               else:
                  total_force += float(line.split()[3])
                  numb_force += 1
#        avg_force[distance] = total_force/float(numb_force)
        if numb_force ==0:
           print distance, "too short simulaiton"
        else:
           avg_force.append(total_force/float(numb_force))
        md_steps.append(numb_force)
        input_force.close()
    
    if len(arg)>3:
       f = open('distance.dat','w')
       f.write("#  %s\n"%('  '.join([distance for distance in distances])))
       for i in range(min([len(atom_distance[distance]) for distance in distances])):
           f.write("%d %s\n"%(i, '   '.join(str(atom_distance[distance][i]) for distance in distances)))
       f.close()
    output.write("# distance   force   free_energy(eV)  md_steps\n")
    for i in range(len(avg_force)):
        if i == 0:
           free_energy.append(avg_force[i])
        else:
           #free_energy.append(simpson(avg_force, float(distance[i]), float(distance[))
           free_energy.append((avg_force[i]+avg_force[i-1])*27.21*(float(distances[i])-float(distances[i-1]))/(2*0.529)+free_energy[i-1]) 
        output.write("%s %12.8f %12.8f %d\n"%(distances[i], avg_force[i], -free_energy[i], md_steps[i]))
if __name__ == '__main__':
    main()
    
