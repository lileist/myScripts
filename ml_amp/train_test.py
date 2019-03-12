from ase.io import read
from ase.io import Trajectory
from amp import Amp
import numpy as np
import sys

arg=sys.argv

configs = Trajectory(arg[1], 'r')

calc = Amp.load('amp.amp')
e_off_log = open('e_off_log.dat','w')
f_off_log = open('f_off_log.dat','w')
e_off_log.write("{:12s} {:12s}\n".format('DFT_E', 'AMP_E'))
f_off_log.write("{:12s} {:12s}\n".format('DFT_f', 'AMP_f'))
off_xyz = open('f_off_traj.xyz','w')
e_off_xyz = open('e_off_traj.xyz','w')

e_log = open('e_log.dat', 'w')
f_log = open('f_log.dat', 'w')
e_log.write("{:12s} {:12s}\n".format('DFT_E', 'AMP_E'))
f_log.write("{:12s} {:12s}\n".format('DFT_f', 'AMP_f'))

traj_amp = Trajectory('amp_'+arg[1],'w')

f_threshold = 2.0
e_threshold = 1.0
index = 0
for config in configs:
  dft_fs = config.get_forces()
  dft_e = config.get_potential_energy()
  config.set_calculator(calc)
  amp_fs = config.get_forces()
  amp_e = config.get_potential_energy()

  traj_amp.write(config)

  e_log.write("{:6d} {:12.6f} {:12.6f} \n".format(index, dft_e, amp_e[0]))
 
  for dft_f, amp_f in zip(dft_fs, amp_fs):
     for i in range(3):
       f_log.write("{:12.6f} {:12.6f} \n".format(dft_f[i], amp_f[i]))

  if np.any(abs(amp_fs - dft_fs) > f_threshold):
     f_off_log.write("{:6d} {:12.6f} {:12.6f} \n".format(index, dft_e, amp_e[0]))
     off_xyz.write("{:4d}\n".format(len(config)))
     off_xyz.write("DFT:{:12.8f} Amp:{:12.8f}\n".format(dft_e, amp_e[0]))
     for i in range(len(config)):
        off_xyz.write("{:4s} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} \n".format(\
                        config[i].symbol, config[i].x, config[i].y, config[i].z, \
                        dft_fs[i][0], dft_fs[i][1], dft_fs[i][2],\
                        amp_fs[i][0], amp_fs[i][1], amp_fs[i][2]))

  if abs(amp_e[0] - dft_e) > e_threshold:
     e_off_log.write("{:6d} {:12.6f} {:12.6f} \n".format(index, dft_e, amp_e[0]))
     e_off_xyz.write("{:4d}\n".format(len(config)))
     e_off_xyz.write("DFT:{:12.8f} Amp:{:12.8f}\n".format(dft_e, amp_e[0]))
     for i in range(len(config)):
        e_off_xyz.write("{:4s} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} \n".format(\
                        config[i].symbol, config[i].x, config[i].y, config[i].z, \
                        dft_fs[i][0], dft_fs[i][1], dft_fs[i][2],\
                        amp_fs[i][0], amp_fs[i][1], amp_fs[i][2]))
  index+=1
  e_log.flush()
  f_log.flush()
  e_off_log.flush()
  f_off_log.flush()
  off_xyz.flush()
  e_off_xyz.flush()
