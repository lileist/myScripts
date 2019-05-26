#!/usr/bin/env python
"""
This code is used to calculated the strain of particels with eon code and lammps pot by linearly interpolating images between transition state and reactant state. The atoms in the reacting center will be reset to their coordinates in the reactant state. 'fix_atoms' is used to define the boundary for strain calculation.
Example input:
max_n = 5
n_image = 8
reacting_center = 1 2 3 4
fix_atoms = 5 6
"""

import sys
import os
from ase.io import read,write
from ase.constraints import FixAtoms
import numpy as np

#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)



def find_area(atoms, index):
    return atoms.get_distance(index[0], index[1]) * atoms.get_distance(index[1], index[2])
#   x_matrix = [atoms[i].x for i in fix_atoms]
#   z_matrix = [atoms[i].z for i in fix_atoms]
#   return (max(x_matrix) - min(x_matrix)) * (max(z_matrix) - min(z_matrix))
   
def readinputs(filename):
    f=open(filename, 'r')
    parameters = {}
    lines=f.readlines()
    for line in lines:
      if line.startswith('#'):
         continue
      fields = line.split('=')
      parameters[fields[0].strip()]=fields[1].replace("\n","").strip()
    return parameters

args = sys.argv
#transition state
paras = readinputs(args[1])
max_n = int(paras['max_n'])
#number of images inerted
n_image = int(paras['n_image'])

workdir = os.getcwd()

configs = []
reacting_center = np.array([int(field) for field in paras['reacting_center'].split()])
fix_atoms =  np.array([int(field) for field in paras['fix_atoms'].split()])


rs = read(str(0)+'.con')
ts = read(str(max_n)+'.con')
n_Pd = len([atom.index for atom in rs if atom.symbol=='Pd'])
n_Au = len(rs) - n_Pd

#rs_area = poly_area(np.array([rs[i].position for i in fix_atoms[0:4:]])) + poly_area(np.array([rs[i].position for i in fix_atoms[2:6:]]))
#rs_area = poly_area(np.array([rs[i].position for i in fix_atoms])) 
rs_area = find_area(rs, fix_atoms)
diff_area = []
#construct dx
disp = (ts.get_positions() - rs.get_positions())/float(n_image)

for i in range(n_image):
   configs.append(rs.copy())

c = FixAtoms(fix_atoms)
for i in range(0, n_image):
  configs[i].set_positions(rs.get_positions()+disp * i)
  for j in reacting_center:
    configs[i][j].position = rs[j].position
  configs[i].set_constraint(c)
#  diff_area.append(find_area(configs[i], fix_atoms)-rs_area)
#  diff_area.append(poly_area(np.array([configs[i][j].position for j in fix_atoms[0:4:]])) + poly_area(np.array([configs[i][j].position for j in fix_atoms[2:6:]])) -rs_area)
#  diff_area.append(poly_area(np.array([configs[i][j].position for j in fix_atoms]))  -rs_area)
  diff_area.append(find_area(configs[i], fix_atoms)-rs_area)
  write(str(i)+'.con',configs[i])

ds = np.linalg.norm(disp[18])
print ds

energy = []
e_out = open(workdir+'/e_strain.dat','w')
for i in range(n_image):
  os.chdir(workdir)
  if os.path.exists(str(i)):
     os.system('rm -rf '+str(i))
  os.makedirs(str(i))
  inputfile = open(str(i)+'.con','r')
  output = open(str(i)+'/pos.con','w')
  n_line = 0
  pd = []
  au = []
  while True:
     line = inputfile.readline()
     n_line += 1
     if not line:
       break
     if n_line == 8 or n_line == 9:
        output.write('%s %s\n'%(line.split()[1], line.split()[0]))
     if n_line < 8:
        output.write(line)
     if 'Pd' in line:
        pd.append(line)
        inputfile.readline()
        pd.append("Coordinates of Component 1\n")
        for j in range(n_Pd):
          pd.append(inputfile.readline())
     if 'Au' in line:
        au.append(line)
        inputfile.readline()
        au.append("Coordinates of Component 2\n")
        for j in range(n_Au):
          au.append(inputfile.readline())
  for line in pd:
     output.write(line)
  for line in au:
     output.write(line)
  output.close()
  os.system('cp /home/leili/exonmobile-project/neb/core_shell/geo_opt/config.ini '+str(i))
  os.system('cp -r /home/leili/exonmobile-project/neb/core_shell/geo_opt/potfiles '+str(i))
  os.chdir(workdir+'/'+str(i))
  os.system('~/code/eon/bin/eon')
  results = open('./output/results.dat','r')
  energy.append(float(results.readlines()[-1].split()[0]))

strain = []
strain_i =  None
for i in range(n_image):
   strain.append(diff_area[i]/rs_area)
   if strain[i] > 0.05 and not strain_i: 
      strain_i = i
curvature = np.polyfit(diff_area[0:strain_i:], energy[0:strain_i:], 2)
curvature_s = np.polyfit(strain[0:strain_i:], energy[0:strain_i:], 2)
#fitting in all range
#curvature = np.polyfit(diff_area, energy, 2)
e_out.write("#f=a*x^2+bx+y to area  : %12.6f %12.6f %12.6f\n"%(curvature[0], curvature[1], curvature[2]))
e_out.write("#f=a*x^2+bx+y to strain: %12.6f %12.6f \n"%(curvature_s[0], curvature_s[0]/rs_area))
for i in range(n_image):
  e_out.write("%8.6f %6.4f %12.6f %12.6f\n"%(diff_area[i], strain[i], energy[i], energy[i]-energy[0]))
