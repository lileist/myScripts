from ase.io import read
import sys
import numpy as np
from ase.io import Trajectory

args = sys.argv
images = read(args[1], index=":") 
opt_traj = Trajectory('center_opt.traj', 'w')
cm_o = None
interval = int(len(images)/ 2000)
for i in range(len(images)):
   atoms = images[i]
   if cm_o is None:
      cm_o = atoms.get_center_of_mass()
   if i % interval == 0:
      cm = atoms.get_center_of_mass()
      atoms.translate(cm_o - cm)
      opt_traj.write(atoms)
