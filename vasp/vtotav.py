#!/usr/bin/env python
"""
This code is used to find NELECT to match electrode potential to specific value relative to SHE (4.43 eV).
Method: calculate Ef (Fermi Energy relative to vaccum level), then U (electrode potential = (-Ef-SHE)/e)
Steepest desent algorithm used to minimize dU (U-U_target)
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
import numpy
import subprocess

sys.stdout.flush()

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def vaccum_pot(idir, index):
#  spin = int(arg[2])
  pot_input = open('LOCPOT', 'r')
  while True:
      line = pot_input.readline()
      if not line:
         break
      for i in range(6):
          line = pot_input.readline()
      ions = [int(field) for field in line.split()]
      print "ions:", ions
      n_ions = sum(ions)
      print "total number of ions:", n_ions
      line = pot_input.readline()
      for i in range(n_ions):
          pot_input.readline()
      pot_input.readline()
      vlocal = []
      vlocal.append(0.0)
#      for i in range(spin):
      fields = [int(field) for field in pot_input.readline().split()]
      ngx = fields[0]
      ngy = fields[1]
      ngz = fields[2]
      print "ngx ngy ngz:",fields 
      nplwv = ngx*ngy*ngz
      #first element is not used
      for j in xrange(int(nplwv/5)+1):
          vlocal.extend([float(field) for field in pot_input.readline().split()])
      print len(vlocal)
      if len(vlocal)>=nplwv+1:
         break
  if idir == 'x':
     nout = ngx
  elif idir == 'y':
     nout = ngy
  else:
     nout = ngz
  scale = 1./float(nplwv/nout)
  vav = []
  vav.append(0.0)
  for i in xrange(nout):
     vav.append(0.0)
  if idir == 'x':
    for i in xrange(1, ngx+1):
        for j in xrange(1, ngy+1):
            for k in xrange(1, ngz+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[i] = vav[i]+ vlocal[ipl] * scale
  elif idir == 'y':
    for j in xrange(1, ngy+1):
        for k in xrange(1, ngz+1):
            for i in xrange(1, ngx+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[j] = vav[j]+ vlocal[ipl] * scale
  elif idir == 'z':
    for k in xrange(1, ngz+1):
        for j in xrange(1, ngy+1):
            for i in xrange(1, ngx+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[k] = vav[k]+ vlocal[ipl] * scale
  output = open('loc_pot.dat','w')
  for i in xrange(1, len(vav)):
      output.write("%d %15.6f\n"%(i, vav[i]))
  tvav = 0.0
  for i in range(index[0], index[1]):
      tvav = tvav+vav[i]
      print vav[i]
  return tvav/float(index[1]-index[0])

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
    doscar = open('DOSCAR', 'r')
    for i in range(6):
        line = doscar.readline()
    fermi = float(line.split()[3])
    doscar.close()
    vaccum_level = vaccum_pot(arg[1], index=[int(arg[2]),int(arg[3])])
    print "Fermi   Vaccum level"
    print fermi, vaccum_level
if __name__ == '__main__':
    main()
    
