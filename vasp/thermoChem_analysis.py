#!/usr/bin/env python
"""
This code is used to summarize theromchemistry data:
vibration: calculated with vasp
data processed with "~/bin/myScripts/data_process/thermochemistry.py"
An example of inputfile (any line starts with '#' will be ignored):
    directories = 
    fix_elements = Au Ag
    stru_filename = CONTCAR
    freq_filename = OUTCAR
    geometry = linear
    symmetrynumber = 1
    temperature = 298
    pressure = 101325
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
import cPickle as pickle
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
    #lines= inp.readlines()
    output = open('free_energy.dat','w')
    energy = None
    output.write("%20s  %s %s %s %s %s %s\n"%("dir", "energy", "e_zpe", "cp_dt", "ts", "g_corr", "free_e"))
    geometry = None
    for dir in dirs: 
       if dir.split("_")[0] == "CO":
          geometry = 'linear'
       elif dir.split("_")[0] =="H":
          geometry = 'monatomic'
       elif dir.split("_")[0] == "COOH":
          geometry = 'nonlinear'
       if geometry is None:
          geometry = 'nonlinear'
       with cd(dir):
            proc = subprocess.Popen("grep 'f  =' OUTCAR", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            vasp_check = proc.communicate()
            if "f  =" not in vasp_check[0]:
                energy = None
                print(dir,"is not done")
            else:
                proc_2 = subprocess.Popen("grep 'energy  without entropy' OUTCAR |tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                vasp_output = proc_2.communicate()
                energy = vasp_output[0].split()[6]

                atoms = read(paras['stru_filename'], index = 0, format='vasp')
                #treat all vibrations
                #del atoms[[atom.index for atom in atoms if atom.symbol in fix_elements]]
                commands = "python ~/myCode/myScripts/data_process/thermochemistry.py " \
                           + paras['format']+" " \
                           + paras['stru_filename']+" " \
                           + paras['freq_filename']+" " \
                           + geometry +" " \
                           + paras['symmetrynumber']+" " \
                           + energy+" " \
                           + paras['temperature']+" " \
                           + paras['pressure'] 
                print(commands)
                proc = subprocess.Popen(commands,
                                        shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                proc.stdin.write(pickle.dumps(atoms))
                proc.stdin.flush()
                thermo_output = proc.communicate()
                lines = thermo_output[0].strip().split("\n")
                cp_dt = 0
                ts = 0
                for line in lines:
                   if "E_ZPE" in line:
                      e_zpe = line.split()[1]
                   #if "C_v -> C_p" in line:
                   #   cp_dt +=float(line.split()[3])
                   #if "Cv_trans" in line or "Cv_rot" in line or "Cv_vib" in line:
                   if "Cv_vib" in line:
                      cp_dt += float(line.split()[2])
                   #if "S_rot" in line or "S_elec" in line or "S_vib" in line:
                   if "S_elec" in line or "S_vib" in line:
                      ts += float(line.split()[3])
                   #if "S (1 atm -> P)" in line:
                   #   ts += float(line.split()[7])
                   if "free_energy" in line:
                      free_e = line.split()[1]
                g_corr = cp_dt - ts + float(e_zpe)
       if energy is None:
          free_e = None
       else:
          free_e = str(float(energy)+g_corr)
       print(dir)
       output.write("%20s  %s %s %s %s %s %s\n"%(dir, energy, e_zpe, str(cp_dt), str(0-ts), str(g_corr), free_e))
if __name__ == '__main__':
    main()
    
