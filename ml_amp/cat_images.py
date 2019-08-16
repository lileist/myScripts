from ase.io import Trajectory
import sys
from set_forces import force_setter
import numpy as np

arg = sys.argv
configs_1=Trajectory(arg[1],'r')
configs_2=Trajectory(arg[2],'r')
configs=Trajectory(arg[3],'w')
n_new =0
removed = 0
for config in configs_1:
   calc = force_setter(energy=config.get_potential_energy(), forces=config.get_forces())
   config.set_cell([[40.,0,0],[0,40.,0],[0,0,40.]],scale_atoms=False)
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
print '1st:', removed, n_new
removed = 0
for config in configs_2:
   calc = force_setter(energy=config.get_potential_energy(), forces=config.get_forces())
   config.set_cell([[40.,0,0],[0,40.,0],[0,0,40.]],scale_atoms=False)
   config.set_pbc((True, True, True))
   config.center()
   config.set_calculator(calc)
   config.get_potential_energy()
   fs=config.get_forces()
   if np.amax(np.absolute(fs)) > 10.:
      removed+=1
      continue
   configs.write(config)
print '2nd:', removed
