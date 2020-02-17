from ase.thermochemistry import HinderedThermo
import numpy
from ase.io import read
#from ase.optimize import QuasiNewton
#from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
import subprocess
import sys
"""
                      1          2               3             4        5                6   7   8
thermochemistry.py [format] [pos_filename] [freq_filename] [nonlinear] [symmetrynumber] [E] [T] [P]
"""

def read_vibs(filename, format = 'vasp'):
    vibs = []
    if format == 'vasp':
       label = 'f  ='
    if format == 'gaussian':
       label = 'Frequencies --'
    cmds = "grep "+"'"+label+"'"+' '+filename
    proc = subprocess.Popen("grep "+"'"+label+"'"+' '+filename, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE) 
    output = proc.communicate()
    for line in output[0].strip().split('\n'):
      if format == 'vasp':
         vibs.append(float(line.split()[7]))
      if format == 'gaussian':
         fields = line.split()
         vibs.append(float(fields[2]))
         vibs.append(float(fields[3]))
         vibs.append(float(fields[4]))
    #print('vibs:', vibs)
    return numpy.array(vibs) 
       
   
#proc_2 = subprocess.Popen("grep 'energy  without entropy' OUTCAR |tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
#vasp_output = proc_2.communicate()
#energy = vasp_output[0].split()[6]

arg = sys.argv
if arg[1] == 'gaussian':
   fmt = 'xyz'
if arg[1] == 'vasp':
   fmt = 'vasp'
atoms = read(arg[2],index= 0,format=fmt)
potentialenergy = float(arg[6])

#print('arg3',arg[3])
vibs = read_vibs(arg[3],arg[1])     #read vibs from file in array format
#print(vibs)
vib_energies = vibs / 8065.54429  # convert to eV from cm^-1
#vib_energies = vib.get_energies()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry=arg[4],
                        symmetrynumber=int(arg[5]), spin=0)
G = thermo.get_gibbs_energy(temperature=float(arg[7]), pressure=float(arg[8]))
