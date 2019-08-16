import sys, random, copy
from ase.io import Trajectory
import numpy as np

arg = sys.argv
configs = Trajectory(arg[1],'r')
traj = Trajectory('train_images.traj', 'w')
test_traj = Trajectory('test_images.traj', 'w')
#del configs[:min_index]

i=0
e_log=[]
configs_list=[]
#f_threshold = 10
f_threshold = 1000
f_split = 100000
f_log =open('f_log.dat', 'w')
f_numb = 0
high_f = []
for config in configs:
  i+=1
  curr_f = config.get_forces()
  #if np.any(abs(curr_f)) > f_threshold:
  if np.amax(np.absolute(curr_f)) > f_threshold:
     continue
  if np.amax(np.absolute(curr_f)) > f_split:
     high_f.append(config)
     continue
  configs_list.append(config)
  e_log.append(config.get_potential_energy())
  for dft_f in curr_f:
     for i in range(3):
       f_numb += 1
       f_log.write("{:10d} {:12.6f} \n".format(f_numb,dft_f[i]))

print len(high_f)
random.shuffle(configs_list)
n_split = int(0.8 * (len(configs_list) + len(high_f))) - len(high_f)
images_train = configs_list[:n_split]
images_test = configs_list[n_split:]
images_train.extend(high_f)

test = open('test_e.dat','w')
train = open('train_e.dat','w')
i=0
for image in images_test:
    i+=1
    test.write("%10d %12.8f\n"%(i, image.get_potential_energy()))
    test_traj.write(image)
i=0
for image in images_train:
    i+=1
    train.write("%10d %12.8f\n"%(i, image.get_potential_energy()))
    traj.write(image)
