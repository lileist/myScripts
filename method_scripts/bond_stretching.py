#!/usr/bin/env python
"""
This code is used to stretch the bond.
strain_cal_cluster [inputfile]
An example of inputfile (any line starts with '#' will be ignored):
    start_structure = POSCAR
    start_atom = 0
    end_atom   = 1
    shift_targets = 2 3
    dstrain    = 0.0025
    njobs      = 10
    fix_direction = all
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

def interpolate(start_atom, end_atom, dstrain, njob=10):
    #shift the start_atom to the original point (end_atom shifts in the same way), then shift end_atom by (i*dstrain) in all direction
    #then shift them back (add back (ax, ay, az))
    dr = []
    (ax, ay, az)=start_atom.position
    (x0, y0, z0)=end_atom.position
    for i in range(njob):
    #    new_coords.append([(x0-ax)*(1 + i * dstrain)+ax, 
    #                       (y0-ay)*(1 + i * dstrain)+ay, 
    #                       (z0-az)*(1 + i * dstrain)+az])
        dr.append([(x0-ax)*(1 + i * dstrain) + ax - x0, 
                   (y0-ay)*(1 + i * dstrain) + ay - y0, 
                   (z0-az)*(1 + i * dstrain) + az - z0])
    #new_coords.append(end_atom.position)
    return dr

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    print paras
    start_atom = int(paras['start_atom'])
    end_atom = int(paras['end_atom'])
    shift_targets = numpy.array([int(field) for field in paras['shift_targets'].split()])
    dstrain = float(paras['dstrain'])
    fixed_direction = paras['fix_direction']
    #read geometry
    p1 = read(filename=paras['start_structure'], index=0, format = 'vasp')
    bond_length = p1.get_distance(start_atom, end_atom)
    fixed_atoms = constrained_indices(p1)
    working_dir = os.getcwd()
    #move atoms based on the given strain and submit jobs
    drs = interpolate(p1[start_atom], p1[end_atom], dstrain, int(paras['njobs']))
    #print new_atoms
    i = 0
    log = open('stretch.dat','w')
    for dr in drs:
        for target in shift_targets:
          p1[target].position += dr
        job_i = working_dir+'/'+str(i)
        make_dir(job_i)   
        write(filename=job_i+'/POSCAR', images=p1, format='vasp')
        fix_atoms(filename=job_i+'/POSCAR', fixed_atoms = [end_atom], fixed_direction=fixed_direction)
        log.write("%4d %8.6f\n"%(i, i*dstrain*bond_length))
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
    
