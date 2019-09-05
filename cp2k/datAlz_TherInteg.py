#!/usr/bin/env python
"""
This code is used to integrate the mean force from thermal integration, thus obtaining free energy
datAlz_TherInteg.py [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
distance = 2.4  2.5  2.6  2.8  3.0  3.2  3.4 3.5 3.6  3.8  4.0  4.2  4.4  4.6  4.8 5.0
mode = accumulative or separate
cp2k_inp = suppl.inp
job_submit_script = run-cp2k.sub
job_submit_cmd    = sbatch
stablize_step = 20000
d_step = 4000
max_step = 400000
reactant = 4.2  5.0
product  = 2.4  2.8
dx = 0.1
index_1 = 14
index_2 = 16
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
import numpy as np
import math
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
              positions.append( np.array([ float(fields[j+1]) for j in range(3) ] ))
       #distance.append(la.norm(positions[0]-positions[1]))
       distance.append(math.sqrt((positions[0][0]-positions[1][0])**2+(positions[0][1]-positions[1][1])**2+(positions[0][2]-positions[1][2])**2))
       if line_numb > max_numb:
          break
    print 'images',line_numb
    f.close()
    return np.array(distance)

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    d_step = int(paras['d_step'])
    if 'mode' in paras.keys():
       mode = paras['mode']
    else:
       mode = 'accumulative'
    stablize_steps = int(paras['stablize_step'])
    distances=paras['distance'].split()
    
    rs = paras['reactant'].split()
    fs = paras['product'].split()
    ts = paras['transition_state'].split()
    rs_min = distances.index(rs[0])
    rs_max = distances.index(rs[1])
    fs_min = distances.index(fs[0])
    fs_max = distances.index(fs[1])
    ts_min = distances.index(ts[0])
    ts_max = distances.index(ts[1])


    output = open(mode+'_freeEnergy.dat','w')
    dx = float(paras['dx'])
    forces=OrderedDict()
    avg_forces = OrderedDict()
    free_energies = OrderedDict()
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
              std_data.write("%s %15.6f %15.6f\n"%(distance, np.average(atom_distance[distance][int(paras['stablize_step']):]),np.std(atom_distance[distance][int(paras['stablize_step']):])))
           else:
              std_data.write("%s %15.6f %15.6f\n"%(distance, np.average(atom_distance[distance][int(paras['stablize_step']):int(paras['max_step'])]),np.std(atom_distance[distance][int(paras['stablize_step']):int(paras['max_step'])])))
        for line in lines:
#            numb_line += 1
            if 'Shake  Lagrangian' in line:
               forces[distance].append(float(line.split()[3]))
        forces[distance]=np.array(forces[distance])
        input_force.close()
        n_forces = len(forces[distance])
        md_steps.append(n_forces)
    print md_steps
    n_slices = int((np.amin(md_steps)-stablize_steps)/d_step)+1
    print n_slices

    if mode == 'accumulative':
       print 'accumulative'
       for key in forces.keys():
           for i in range(n_slices):
             end_index = stablize_steps+(i+1)*d_step+1
             if end_index > len(forces[key]):
                end_index = -1
             if key not in avg_forces: 
                #avg_forces[key] = [np.average(forces[key][stablize_steps:stablize_steps+(i+1)*d_step+1])]
                avg_forces[key] = [np.average(forces[key][stablize_steps:end_index])]
             else:
             #avg_forces[key].append(np.sum(forces[key][stablize_steps+i*d_step:stablize_steps+(i+1)*d_step+1)])
                #avg_forces[key].append(np.average(forces[key][stablize_steps:stablize_steps+(i+1)*d_step+1]))
                avg_forces[key].append(np.average(forces[key][stablize_steps:end_index]))
    elif mode == 'separate':
       for key in forces.keys():
           for i in range(n_slices):
             end_index = stablize_steps+(i+1)*d_step+1
             if end_index > len(forces[key]):
                end_index = -1
             if key not in avg_forces: 
                avg_forces[key] = [np.average(forces[key][stablize_steps+i*d_step:end_index])]
             else:
                avg_forces[key].append(np.average(forces[key][stablize_steps+i*d_step:end_index]))
    elif mode == 'step_forward':
       print mode
       for key in forces.keys():
           for i in range(n_slices):
             if key not in avg_forces: 
                avg_forces[key] = [np.average(forces[key][i*d_step:i*d_step+stablize_steps])]
             else:
                avg_forces[key].append(np.average(forces[key][i*d_step:i*d_step+stablize_steps]))
    if len(arg)>3:
       f = open('distance.dat','w')
       f.write("#  %s\n"%('  '.join([distance for distance in distances])))
       for i in range(min([len(atom_distance[distance]) for distance in distances])):
           f.write("%d %s\n"%(i, '   '.join(str(atom_distance[distance][i]) for distance in distances)))
       f.close()
    output_head="{:8s}".format("# CV")
    
    j_key = 0
    """
    free_energies: {'0': a list of free_energy with time from stablize_step to 0*d_step
                     ...
                    'n_slices-1' : a list of free_energy with time from stablize_step to n_slices-1 * d_step
                     }
    """
    fe_slices={}
    reactions=open(mode+'_reaction_es.dat','w')
    reactions.write('{:12s} {:12s} {:12s} {:12s} {:12s} {:12s}\n'.format(
                    'time','rs','ts','fs','reactionE','barrier'))
    for i in range(n_slices):
       fe_slices[i] = []
       for j in range(len(distances)):
           key_1 = distances[j]
           if key_1 not in free_energies:
              free_energies[key_1]=[]
           print key_1
           if j == 0:
              temp_fe=avg_forces[key_1][i]
           else:
              key_0 = distances[j-1]
              #free_energy.append(simpson(avg_force, float(distance[i]), float(distance[))
              temp_fe = (avg_forces[key_1][i]+avg_forces[key_0][i])* \
                          27.21*(float(key_1)-float(key_0)) \
                          /(2*0.529)-free_energies[key_0][i] 
           free_energies[key_1].append(-temp_fe)
           fe_slices[i].append(-temp_fe)
       es=np.array(fe_slices[i])
       rs_e = np.amin(es[rs_min:rs_max+1])
       ts_e = np.amax(es[ts_min:ts_max+1])
       fs_e = np.amin(es[fs_min:fs_max+1])
       if mode=='step_forward':
          reactions.write('{:12.2f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n'.format(
                            float(i*d_step)/2000.0, rs_e, ts_e, fs_e, fs_e - rs_e, ts_e - rs_e))
          continue
       end_index =stablize_steps+(i+1)*d_step
       if end_index > np.amin(md_steps):
          end_index = np.amin(md_steps)
       reactions.write('{:12.2f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n'.format(
                         float(end_index)/2000.0, rs_e, ts_e, fs_e, fs_e - rs_e, ts_e - rs_e))
    for i in range(len(free_energies[distances[0]])):
       output_head += "{:13s}".format('FreeEnergy'+str(i))
    output_head += "{:13s}".format('Average')
    output_head += "{:13s}".format('STD')
    output.write("%s \n"%(output_head))
    for i in range(len(distances)):
       outdata = format_joint(free_energies[distances[i]])
       outdata += "{:12.6f}".format(np.average(free_energies[distances[i]]))
       outdata += "{:12.6f}".format(np.std(free_energies[distances[i]], ddof=1))
       output.write("%8s %s \n"%(distances[i], outdata))

if __name__ == '__main__':
    main()
    
