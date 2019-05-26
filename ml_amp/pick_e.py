import sys, random, copy
from ase.io import Trajectory
import numpy as np

def log_atoms(flag, configs, traj, test_traj):
   images_test = []
   images_train = []
   for i in range(len(flag)):
     if i == len(flag)-1 and len(configs) >flag[-1]:
        temp=configs[flag[-1]:len(configs)]
     else:   
        temp=configs[flag[i]:flag[i+1]]
     try:
        index = len(temp)/2
        selected = temp.pop(index)
        traj.write(selected)
        images_train.append(selected)
        selected = random.choice(temp)
        test_traj.write(selected)
        images_test.append(selected)
     except:
        continue
   return images_test, images_train

arg = sys.argv
configs = Trajectory(arg[1],'r')
min_index=int(arg[2])
traj = Trajectory("train_"+arg[1].split('.')[0]+'.traj', 'w')
test_traj = Trajectory("test_"+arg[1].split('.')[0]+'.traj', 'w')
min_e = configs[min_index].get_potential_energy()
#del configs[:min_index]

split_1= -56.70

#min_e-->split_1-->max_e
low_split_flag=[]
high_end=[]
n_low = 3000
#n_high = 
#create split criteria for low_end
de_low   = (split_1-min_e)/float(n_low)
for j in range(n_low):
   #print selected
   curr_e = min_e + j*de_low
   if curr_e >split_1:
      break
   low_split_flag.append(curr_e)

i=0
e_log=[]
configs_list=[]
for config in configs:
  i+=1
  curr_e = config.get_potential_energy()
  if i<min_index+1 or curr_e > -40.0:
    continue
  if curr_e > split_1:
     high_end.append(config)
     continue
  configs_list.append(config)
  e_log.append(config.get_potential_energy())

flag=np.searchsorted(e_log, low_split_flag)
images_test, images_train = log_atoms(flag, configs_list, traj,test_traj)
print "high end:", len(high_end)
for config in high_end:
  traj.write(config)
  images_train.append(config)

test = open('test_e.dat','w')
train = open('train_e.dat','w')
i=0
for image in images_test:
    i+=1
    test.write("%10d %12.8f\n"%(i, image.get_potential_energy()))
i=0
for image in images_train:
    i+=1
    train.write("%10d %12.8f\n"%(i, image.get_potential_energy()))
