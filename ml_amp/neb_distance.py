from ase.io import read
import numpy as np

images = read('movie.traj',index=":",format='traj')
rs= images[0].get_positions()
for i in range(len(images)):
   print np.linalg.norm(images[i].get_positions() - rs)
