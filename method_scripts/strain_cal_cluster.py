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
import numpy
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
    new_atoms=[]
    new_atoms.append(start_atom)
    z = start_atom.z
    x = start_atom.x
    for i in range(n-1):
        x_new = x+(n+1)*dx
        new_atoms.append(start_atom.position=(x_new, k*x_new+b, z))
    return new_atoms.append(end_atom)

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    print paras
    boundary_atoms = numpy.array([float(field) for field in paras['boundary_atoms'].split()])
    #numpy.array([0, 1, 2, 38, 39])
    boundary = numpy.array([float(field) for field in paras['boundary'].split()])
    #([0, 4.407, 0])
    strain = numpy.array([float(field) for field in paras['strain'].split()])
    strain_direction = paras['strain_direction']
    span = float(paras['span'])
    active_space = ([field for field in paras['active_space'].split()])
    #initialize strain for each direction
    dx = numpy.zeros(len(strain))
    dy = dx
    dz = dx
    #set strain
    if strain_direction == 'x':
       dx = strain
       fixed_direction = 1
    if strain_direction == 'y':
       dy = strain
       fixed_direction = 2
    if strain_direction == 'z':
       dz = strain
       fixed_direction = 3
    #read geometry
    p1 = read(filename=paras['structure_file'], index=0, format = 'vasp')
    fixed_atoms = constrained_indices(p1)
    working_dir = os.getcwd()
    #move atoms based on the given strain and submit jobs
    for i in range (0,len(strain)):
        p2 = p1.copy()
        for atom in p2:
            if atom.index in fixed_atoms or atom.symbol not in active_space:
               continue
            #move atoms at the boundary. Keep the original shape of the boundary
            if atom.index in boundary_atoms:
               atom.position = atom.position +(span*dx[i],span*dy[i], span*dz[i])
               continue
            #move other atoms. 
            atom.position = tuple(boundary + (numpy.array(atom.position) - boundary)*numpy.array([1+dx[i],1+dy[i],1+dz[i]]))
        if strain[i]<0:
           job_i = working_dir+'/'+'c_'+str(abs(strain[i]))
        else:
           job_i = working_dir+'/'+str(strain[i])
        make_dir(job_i)   
        write(filename=job_i+'/POSCAR', images=p2, format='vasp')
        fix_atoms(filename=job_i+'/POSCAR', fixed_atoms = boundary_atoms, fixed_direction=fixed_direction)
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
    
