import sys
import numpy as np

args = sys.argv
fin = open(args[1], 'r')
fp_type = args[2]
if fp_type == 'G2':
  etas = ' '
  rss = ' '
  for line in fin.readlines():
     fields = [field for field in line.split()]
     etas  += "{:8.2f}".format(float(fields[3])) +' '
     rss += "{:8.2f}".format(float(fields[4]))+ ' '
  fout = open('input.ini', 'w')
  fout.write('Rc = 6.0\n')
  fout.write('etas = '+etas+'\n')
  fout.write('Rs = '+rss+'\n')
if fp_type == 'G4':
  zetas = ' '
  thetas = ' '
  gammas = ' '
  for line in fin.readlines():
     fields = [field for field in line.split()]
     zetas  += "{:8.2f}".format(float(fields[3])) +' '
     gamma   = float(fields[5])
     gammas += "{:8.2f}".format(gamma)+ ' '
     #if gamma >0 :
     thetas += "{:8.2f}".format(float(fields[4]))+ ' '
     #if gamma <0 :
     #  thetas += "{:8.2f}".format(-np.pi+float(fields[4]))+ ' '
  fout = open('input.ini', 'w')
  fout.write('Rc = 6.0\n')
  fout.write('zetas = '+zetas+'\n')
  fout.write('theta_s = '+thetas+'\n')
  fout.write('gammas = '+gammas+'\n')
  fout.write('eta = '+'0.005'+'\n')
