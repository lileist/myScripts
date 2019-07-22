#!/usr/bin/env python
"""
This code is used to integrate the mean force from thermal integration, thus obtaining free energy
datAlz_TherInteg.py [inputfile]
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
from collections import OrderedDict 
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

def format_joint(elements):
    jointed = ''
    for element in elements:
       jointed+="{:12.6f}".format(element)
    return jointed

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
    d_step = float(paras['d_step'])
    stablize_steps = int(paras['stablize_step'])
    distances=paras['distance'].split()
    output = open('freeEnergy.dat','w')
    dx = float(paras['dx'])
    forces=OrderedDict()
    avg_forces = OrderedDict()
    free_energies = []
    md_steps = []
    atom_distance = {}
    std_data = open('avg_std.dat', 'w')
    force_file = None
    for distance in distances:
        forces[distance]=[]
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
               forces[distance].append(float(line.split()[3]))
        np.arraye(forces[distance])
        input_force.close()
        n_forces = len(forces[distance])
        md_steps.append(n_forces)

    n_slices = int((np.amin(md_steps)-stablize_steps)/d_step)
    for key in forces.keys():
        for i in range(n_slices):
          #avg_forces[key].append(np.sum(forces[key][stablize_steps+i*d_step:stablize_steps+(i+1)*d_step+1)])
          avg_forces[key].append(np.average(forces[key][stablize_steps:stablize_steps+(i+1)*d_step+1)])
    
    if len(arg)>3:
       f = open('distance.dat','w')
       f.write("#  %s\n"%('  '.join([distance for distance in distances])))
       for i in range(min([len(atom_distance[distance]) for distance in distances])):
           f.write("%d %s\n"%(i, '   '.join(str(atom_distance[distance][i]) for distance in distances)))
       f.close()
    output.write("# distance   force   free_energy(eV)  md_steps\n")
    
    j_key = 0
    outdata=OrderedDict()
    sum_slices=OrderedDict()
    """
    free_energies: {'0': a list of free_energy with time from stablize_step to 0*d_step
                     ...
                    'n_slices-1' : a list of free_energy with time from stablize_step to n_slices-1 * d_step
                     }
    """
    for i in range(n_slices):
       for j in range(len(distances)):
           key_1 = distances[j]
           if j == 0:
              outdata[j]="{:12.6f}".format(avg_forces[key_1][i])
              sum_slices[j] = avg_forces[key_1][i]
              free_energies[key_1].append(avg_forces[key_1][i])
           else:
              key_0 = distances[j-1]
              #free_energy.append(simpson(avg_force, float(distance[i]), float(distance[))
              free_energies[key_1].append((avg_forces[key_1][i]+avg_forces[key_0][i])* \
                                         27.21*(float(key_1)-float(key_0)) \
                                         /(2*0.529)+free_energies[key_0][i]) 
              outdata[j] += "{:12.6f}".format(free_energies[i][j]) 
              sum_slices[j] += free_energies[i][j] 

    for i in range(len(distances)):
       outdat = format_joint(free_energies[distances[i]])
       outdata += "{:12.6f}".format(np.average(free_energies[distances[i]]))
       outdata += "{:12.6f}".format(np.std(free_energies[distances[i]], ddof=1))
       output.write("%s %s \n"%(distances[i], outdata)

if __name__ == '__main__':
    main()
    
