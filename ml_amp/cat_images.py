"""
This code is used to cat a series of trajectory files
"""

from ase.io import Trajectory
import sys
from set_forces import force_setter
import numpy as np
import argparse

arg = sys.argv

parser = argparse.ArgumentParser()
parser.add_argument('--Nfile', type=int, nargs='+', metavar='8', 
                    default=None,
                    help='number of trajectory files')
parser.add_argument('--prefix', type=str, nargs='+', metavar='f', 
                    default='f',
                    help='prefix of trajectory files')
parser.add_argument('--output', type=str, nargs='+', metavar='train.traj', 
                    default='train.traj',
                    help='number of trajectory files')
args = parser.parse_args()
configs=Trajectory(args.output[0],'w')
prefix = args.prefix[0]
for i in range(args.Nfile[0]):
  configs_1=Trajectory(prefix+str(i)+'.traj','r')
  n_new =0
  removed = 0
  for config in configs_1:
     calc = force_setter(energy=config.get_potential_energy(), forces=config.get_forces())
     config.set_cell([[80.,0,0],[0,80.,0],[0,0,80.]],scale_atoms=False)
     config.set_pbc((True, True, True))
     config.center()
     config.set_calculator(calc)
     e=config.get_potential_energy()
     fs=config.get_forces()
     if np.amax(np.absolute(fs)) > 10.:
        removed+=1
        continue
     configs.write(config)
     n_new+=1
  print("{:2d}: {:10d} {:10d}".format(i, removed, n_new))
