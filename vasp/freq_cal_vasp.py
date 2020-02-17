#!/usr/bin/env python
"""
This code is used to find NELECT to match electrode potential to specific value relative to SHE (4.43 eV).
Method: calculate Ef (Fermi Energy relative to vaccum level), then U (electrode potential = (-Ef-SHE)/e)
Steepest desent algorithm used to minimize dU (U-U_target)
An example of inputfile (any line starts with '#' will be ignored):
    directories = 
    fix_elements = Au Ag
    submit = 1
    job_submit_script = qsub.hf
    job_submit_cmd    = sbatch
"""
import sys
import os
import errno
import ase
from ase import Atoms
from ase.io import read, write
import numpy
import subprocess
from ase.constraints import FixAtoms
sys.stdout.flush()

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):

        print(newPath)
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

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    #inp = open(paras['incar'], 'r')
    dirs = paras['directories'].split()
    print(dirs)
    fix_elements = paras['fix_elements']
    try:
      fix_indices = [int(field) for field in paras['fix_indices'].split()]
    except:
      fix_indices =None
      pass
    #lines= inp.readlines()
    output = open('energy.dat','a')
    energy = None
    for dir in dirs:
       with cd(dir):
            proc = subprocess.Popen("grep 'reached required accuracy' OUTCAR", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            vasp_check = proc.communicate()
            print(vasp_check)
            if b"stopping structural energy minimisation" not in vasp_check[0]:
                energy = None
                print(dir,"is not done")
            else:
                proc_2 = subprocess.Popen("grep 'energy  without entropy' OUTCAR |tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                vasp_output = proc_2.communicate()
                energy = vasp_output[0].split()[6]
                os.system('cp OUTCAR OUTCAR_ener')
                p1 = read('CONTCAR', index = 0, format='vasp')
                
                os.system('cp ../INCAR_freq INCAR')
                indices=[atom.index for atom in p1 if atom.symbol in fix_elements]
                if fix_indices is not None:
                   indices.extend(fix_indices)
                print(indices)
                c = FixAtoms(indices)
                p1.set_constraint(c)
                write('POSCAR', p1, format = 'vasp')
            if int(paras['submit']) == 1 and energy is not None:
               os.system(paras['job_submit_cmd']+' '+paras['job_submit_script'])
               print(os.getcwd(),'is submitted')
       output.write("%20s    %s\n"%(dir, energy))
if __name__ == '__main__':
    main()
    
