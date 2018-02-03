#!/usr/bin/env python
import sys
sys.path.append("../../")
import os

#from ase.test import NotAvailable
from ase.io import read
from ase.atoms import Atoms
import numpy
from expectra.cal_exafs import Expectra
#from expectra.basin import make_dir
from expectra.cal_exafs import match_x
from expectra.fft import hanning_window
#prefix 'old' used to destinguish configs in databased and new one
def calc_area(y_exp, y_theory, calc_type='area', average = False):
    if len(y_exp) != len(y_theory):
        print "Warning: number of points in chi_exp and chi_theory is not equal"
        
    numb = min(len(y_exp),len(y_theory))
    area_diff = 0.00
    for i in range(0, numb):
      diff = numpy.absolute(y_exp[i] - y_theory[i])
      area_diff = area_diff + diff
    #print ('%s: %15.6f' % ("area_diff", area_diff))
    return area_diff

def output_exafs(filename,k,chi):
    exafs_out = open(filename, 'w')
    for i in range(len(k)):
       exafs_out.write("%12.6f %12.6f\n"%(k[i],chi[i]))
    exafs_out.close()

exafs_filename = 'all_exafs.dat'
f_exafs = open(exafs_filename,'r')
kmin=2.5
kmax=11.0
kweight = 2
dk = 0.05
rmin = 1
rmax = 3
ft_part = 'mag'

output = open('r_diff.dat','w')

exp_r = open('exp_r.dat','r')
r_exp = []
mag_exp = []
lines = exp_r.readlines()
for line in lines:
    fields=line.split()
    r_exp.append(float(fields[0]))
    mag_exp.append(float(fields[1]))
exp_r.close() 
r_exp, mag_exp = match_x(r_exp,mag_exp,r_exp,rmin,rmax)
output_exafs('scaled_exp_r.dat', r_exp, mag_exp)

exp_chi = open('AuPt_1ML_Au_k2.dat','r')
lines = exp_chi.readlines()
k_exp = []
chi_exp = []
for line in lines:
    fields=line.split()
    k_exp.append(float(fields[0]))
    chi_exp.append(float(fields[1]))
k_exp = numpy.array(k_exp)
chi_exp = numpy.array(chi_exp)

chi_exp /= k_exp**2
window = hanning_window(k_exp, kmin, kmax, dk)
chi_exp *= window
chi_exp *= k_exp**2
k_exp, chi_exp = match_x(k_exp,chi_exp,k_exp,kmin,kmax)
output_exafs('scaled_exp.dat', k_exp, chi_exp)

chi_differ = open('chi_diff.dat','r')

d_chi=[]
lines = chi_differ.readlines()
for line in lines:
   d_chi.append(float(line.split("chi_differ:")[1]))
chi_differ.close()


image_n = 0
while True:
    k=[]
    chi=[]
    image_n += 1

    line = f_exafs.readline()
    if not line or image_n >100:
       break
    exafs_out = open("chi.dat",'w')
    for i in range(382):
        line = f_exafs.readline()
        #print line
        fields=line.split()
        k.append(float(fields[0]))
        chi.append(float(fields[1]))
    k = numpy.array(k)
    chi = numpy.array(chi)

    for i in range(len(k)):
       #print k[i],chi[i]
       exafs_out.write("%12.6f %12.6f\n"%(k[i],chi[i]))
    exafs_out.close()

    window = hanning_window(k, kmin, kmax, dk)
    print window
    chi *= window
    chi *= k**2
    k_scaled, chi_scaled = match_x(k_exp,chi,k,kmin,kmax)
    output_exafs('scaled_theory.dat', k_scaled, chi_scaled)
    area_k = calc_area(chi_exp, chi_scaled)

    xafsft_para = ['xafsft',
                   '--kmin', str(kmin),
                   '--kmax', str(kmax),
                   '--kweight', str(kweight),
                   '--dk', str(dk),
                   '--rmin', str(rmin),
                   '--rmax', str(rmax),
                   '--ft-part', ft_part,
                   'chi.dat']
    join_symbol = ' '
    xafsft_cmd = join_symbol.join(xafsft_para)
    print image_n,": ", xafsft_cmd
    os.system(xafsft_cmd)

    inputfile = 'exafs.chir'
    exafs_r = open(inputfile, 'r')
    lines = exafs_r.readlines()
    r=[]
    mag=[]
    for line in lines:
        fields=line.split()
        r.append(float(fields[0]))
        mag.append(float(fields[1]))
    exafs_r.close()

    r_scaled, mag_scaled = match_x(r_exp,mag,r,rmin,rmax)
    output_exafs('scaled_r.dat', r_scaled, mag_scaled)
    area = calc_area(mag_exp, mag_scaled)
    output.write("%6d %12.6f %12.6f %12.6f\n"%(image_n, d_chi[image_n],area_k,area))
f_exafs.close()

