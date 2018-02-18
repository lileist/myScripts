#!/usr/bin/env python
import sys
import numpy as np
from ase.io import read,write
from ase.build import add_adsorbate
from math import pi

"""
Example of input file:
===================================
slab_file        =  POSCAR             # initial structure
format           =  vasp                  # file format
adsorbate_file   = POSCAR_COOH
#rot_angle       =  -60              # unit in degree
rot_axis         =  z
rot_center       = 0
rotate           = True
non_linear       = True
bidentate        = True
===================================
1.5   180    #height  rot_angle '-' clockwise
1.4,7.5   4.3,7.5
===================================
"""
def readinputs(filename):
    f=open(filename, 'r')
    parameters = {}
    lines=f.readlines()
    for line in lines:
      if line.startswith('#') or len(line.split()) == 0:
         continue
      fields = line.partition('#')[0].split('=')
      if fields[1].replace("\n","").strip() == 'True':
         parameters[fields[0].strip()] = True
         continue
      if fields[1].replace("\n","").strip() == 'False':
         parameters[fields[0].strip()] = False
         continue
      parameters[fields[0].strip()]=fields[1].replace("\n","").strip()
    return parameters

def get_height_sites(filename):
    """
    height_sites =[[height,angle,[index_1,index_2]],[...]]
    """
    f=open(filename, 'r')
    height_sites=[]  
    while True:
      height_site=[]  
      line = f.readline()
      if not line:
         break
      if line.startswith('#'):
         continue
      height_site.append(float(line.split()[0]))
      height_site.append(float(line.split()[1]))
      line = f.readline()
      #fields = [field.replace("\n","") for field in line.split()]
      #print fields
      #for field in fields:
      #    height_site.append(tuple([float(temp) for temp in field.split(',')]))
      height_site.append([int(temp) for temp in line.split()])
      height_sites.append(height_site)
    print height_sites
    return height_sites

#def rotate(adsorbate,ads_vect, vect):
#    rot_axis = numpy.cross(ads_vect,vect)
#    angle = np.linalg.norm(rot_axis)/(np.linalg.norm(ads_vect)*np.linalg.norm(vect))

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    p = read(paras['slab_file'],index=0,format=paras['format'])
    temp = p.copy()
    adsorbate=read(paras['adsorbate_file'],index=0,format='vasp')
    adsorbate.set_cell(p.get_cell())
    temp_ads = adsorbate.copy()
    #for latest ASE version
    #adsorbate.euler_rotate(phi=60.0,   #rotation angle around z
    #                       theta=0.0,  #rotation angle around x
    #                       psi=0.0)    #2nd rotation angle around z
    height_sites = get_height_sites(arg[2])
    if paras['rot_center']:
       center=tuple(adsorbate[int(paras['rot_center'])].position)
       print center
   
    #height_sites = [[1.5,(1.4,7.5),(4.3,7.5)]]
    
    non_linear = paras['non_linear']
    bidentate = paras['bidentate']
    
    i=0
    for height_site in height_sites:
       i+=1
       height = height_site[0]
       indices = height_site[2]
       if height_site[1]!=0:
          adsorbate.rotate(v=paras['rot_axis'],a=height_site[1]*pi/180,center=center)
       if len(indices)==2:
          site = (np.array(p[indices[0]].position) + np.array(p[indices[1]].position))/2
          if non_linear and bidentate:
             if abs(p[indices[1]].z - p[indices[0]].z) > 0.2:
                adsorbate.rotate(v='y',a=-30.0*pi/180,center=center)  #for fcc211
             site = (np.array(p[indices[0]].position) + np.array(site))/2
   
       elif len(indices)==3:
          site = (np.array(p[indices[0]].position) + np.array(p[indices[1]].position) + np.array(p[indices[2]].position))/3
       else:
          site = np.array(p[indices[0]].position)
       print site
       add_adsorbate(p,adsorbate,height,tuple([site[0],site[1]]))
       write('POSCAR'+str(i),images=p,format='vasp')
       p = temp.copy()
       adsorbate = temp_ads.copy()
if __name__ == '__main__':
    main()
