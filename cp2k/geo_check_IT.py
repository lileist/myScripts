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
from ase.io import read, write, Trajectory
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
    tag = None
    atoms =None
    traj = Trajectory('reaction.traj','w', atoms)
    #traj=[]
    cellsize = 20.0
    for distance in distances:
        cp2k_inp = open(distance+'/'+paras['cp2k_inp'], 'r')
        lines = cp2k_inp.readlines()
        for line in lines:
            if 'PROJECT_NAME' in line:
               tag = line.split()[1]
               break
        try:
          p=read(distance+'/'+'geometry.xyz', index=slice(-1, None),format='xyz')[0]
          p.set_cell([[cellsize,0,0],[0,cellsize,0],[0,0,cellsize]],scale_atoms=False)
          p.set_pbc((True, True, True))
        except:
          print(distance)
          continue
        p.center()
        print(p)
        
        #traj.append(p)
        traj.write(p)
        #write('reaction.traj',traj, format='traj')
if __name__ == '__main__':
    main()
    
