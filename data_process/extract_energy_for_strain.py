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

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    print paras
    target_atom = int(paras['target_atom'])
    #boundary_atoms = numpy.array([float(field) for field in paras['boundary_atoms'].split()])
    #numpy.array([0, 1, 2, 38, 39])
    #boundary = numpy.array([float(field) for field in paras['boundary'].split()])
    #([0, 4.407, 0])
    #strain = numpy.array([float(field) for field in paras['strain'].split()])
    image_numb = int(paras['image_numb'])
    #strain_direction = paras['strain_direction']
    #span = float(paras['span'])
    #active_space = ([field for field in paras['active_space'].split()])
    #initialize strain for each direction
    #dx = numpy.zeros(len(strain))
    #dy = dx
    #dz = dx
    """
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
    """
    #read geometry
    #file_pre = os.getcwd().split('/')[-1]
    output = open('stain_e'+'.dat', 'w')
    #move atoms based on the given strain and submit jobs
    int_s = read('./0/CONTCAR', index=0)
    final_s = read('./6/CONTCAR', index=0)
    dist = numpy.linalg.norm(numpy.array(int_s[target_atom].position - final_s[target_atom].position))
    print dist
    dr = dist/float(image_numb+1)
    for i in range (0,image_numb+2):
        #if strain[i]<0:
        #   job_i = './'+'c_'+str(abs(strain[i]))
        #else:
        #   job_i = './'+str(strain[i])
        job_i = str(i)
        grep_energy = "grep 'energy  without entropy' "+job_i + "/OUTCAR" +" |tail -n 1"
        print grep_energy
        #check_completion = "grep 'reached required accuracy' "+job_i + "/OUTCAR" +" |tail -n 1"
        fetch_force = "grep FORCES: "+job_i+"/OUTCAR"+" |tail -n 1"
        p =subprocess.Popen(fetch_force, shell=True, stdout=subprocess.PIPE)
        results = p.communicate()
        print results[0]
        force = results[0].split('RMS')[1].split()[0]
        #if 'reached' in results[0].split()[0]:
        p =subprocess.Popen(grep_energy, shell=True, stdout=subprocess.PIPE)
        results = p.communicate()
        output.write("%6.5f %s %s\n"%(float(i)*dr,results[0].split()[6],force))
        #continue
        #else:
        #   output.write("%6.5f %s\n"%(strain[i],'job not finished'))

if __name__ == '__main__':
    main()
    
