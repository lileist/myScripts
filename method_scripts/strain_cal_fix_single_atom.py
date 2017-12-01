#!/usr/bin/env python
"""
This code is used to calculate strain of nanoparticles.
Method: stretch partilces along certain directions by dx and fix the boundary atoms, do optimization
strain_cal_cluster [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    boundary_atoms= 0 1 2 38 39
    boundary = 0 4.407 0
    strain = -0.03 -0.02 -0.015 -0.01 -0.005 0 0.005 0.01 0.015 0.02 0.03
    strain_direction = y
    span = 7.247
    active_space = Au
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
from ase.constraints import constrained_indices, FixAtoms
import numpy, copy
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

def fix_atoms(filename, fixed_atoms=None, fixed_direction=2):
    posfile = open(filename, 'r')
    lines= posfile.readlines()
    posfile.close()
    output = open(filename, 'w')
    i=0
    for line in lines:
       if i-8 in fixed_atoms:
          fields = line.split()
          if fixed_direction =='all':
             fields[3] = 'F'
             fields[4] = 'F'
             fields[5] = 'F'
          else:
             fields[2+fixed_direction] = 'F'
          line=' '+'  '.join(fields)+'\n'
       output.write("%s"%(line))
       i+=1
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

def interpolate(start_atom, end_atom, n):
    #find slop and intercept
    k = (end_atom.y - start_atom.y)/(end_atom.x - start_atom.x)
    b = end_atom.y - k * end_atom.x
    dx = (end_atom.x - start_atom.x)/float(n)
    print k, b, dx
    new_coords=[]
    new_coords.append(start_atom.position)
    z = start_atom.z
    x = start_atom.x
    dz = (end_atom.z - z)/float(n)
    for i in range(n):
        x_new = x + float(i + 1) * dx
        new_coords.append([x_new, k*x_new+b, z+(i+1)*dz])
    #new_coords.append(end_atom.position)
    return new_coords

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    print paras
    target_atom = int(paras['target_atom'])
    #set strain
    fixed_direction = paras['fix_direction']
    #read geometry
    p1 = read(filename=paras['start_structure'], index=0, format = 'vasp')
    p2 = read(filename=paras['end_structure'], index=0, format = 'vasp')
    fixed_atoms = constrained_indices(p1)
    working_dir = os.getcwd()
    #move atoms based on the given strain and submit jobs
    print p1[target_atom]
    new_atoms = interpolate(p1[target_atom], p2[target_atom], int(paras['image_numb']))
    #print new_atoms
    i = 0
    for atom in new_atoms:
        p1[target_atom].position = atom
        job_i = working_dir+'/'+str(i)
        make_dir(job_i)   
        write(filename=job_i+'/POSCAR', images=p1, format='vasp')
        fix_atoms(filename=job_i+'/POSCAR', fixed_atoms = [target_atom], fixed_direction=fixed_direction)
        i += 1
        os.system('cp POTCAR '+job_i)
        os.system('cp KPOINTS '+job_i)
        os.system('cp INCAR '+job_i)
        os.system('cp '+paras['job_submit_script']+' '+job_i)
        with cd(job_i):
             os.system(paras['job_submit_cmd']+' '+paras['job_submit_script'])
             print os.getcwd(),'is submitted'
        print 'Back to directory',os.getcwd()
if __name__ == '__main__':
    main()
    
