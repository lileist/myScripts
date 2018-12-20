"""
This code is used to read geometries with corresponding energy and force from OUTCAR file
"""
from ase import Atoms
import sys, os, glob
import subprocess
import copy
import numpy as np
from ase.io import Trajectory, read
from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError
from expectra.io import read_con
from expectra.atoms_operator import rot_match


class force_setter(Calculator):
    implemented_properties = ['energy', 'forces']
    default_parameters = {'epsilon': 1.0,
                          'sigma': 1.0,
                          'rc': None}
    nolabel = True

    def __init__(self, energy=None, forces=None,**kwargs):
        Calculator.__init__(self, **kwargs)
        self.energy = energy
        self.forces = forces

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

#        natoms = len(self.atoms)

#        positions = self.atoms.positions
#        cell = self.atoms.cell

        self.results['energy'] = self.energy
        self.results['forces'] = self.forces

def call_subprocess(cmd, shell=True):
    proc = subprocess.Popen(cmd, shell=shell, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
    #proc.stdin.write(pickle.dumps(inputs))
    #proc.stdin.flush()
    output = proc.communicate()
    if output[1] is not None:
       print ("Error {} raised for expectra: {}".format(return_value,output[1]))
       sys.exit()
    return output[0].split('\n')

def read_atoms(dir, elements, shell=True):
    #read dimer traj from climb.con
    dimer_configs = read_con('climb.con')
    if type(dimer_configs) is not list:
       dimer_configs=[dimer_configs]
    #get number of images generated during rs optimization
    cmd="sed -n '/Minimizing initial structure/,/Saddle point search started/p' ll_out |wc -l"
    n_mini_rs = int(call_subprocess(cmd)[0])-2
    #check if dimer is successful
    cmd='Starting Minimization 1'
    n_mini = len(call_subprocess(cmd))
    n_mini_2 = 0
    #if dimer succeed, get number of images generated during min1 and min2 optimization
    if n_mini !=0:
       cmd="sed -n '/Starting Minimization 1/,/Saddle Search/p' ll_out"
       n_mini = len([line for line in call_subprocess(cmd) if 'Matter' in line])
       cmd="sed -n '/Starting Minimization 1/,/Starting Minimization 2/p' ll_out"
       n_mini_1 = len([line for line in call_subprocess(cmd) if 'Matter' in line])
       #if path successfully connected rs and fs, get number of images generated during min2 optimization
       #merge min1 trajectory to rs trajectory
       cmd="grep 'Final status: Success' ll_out"
       if 'Success' in call_subprocess(cmd)[0]:
          n_mini_2 = n_mini - n_mini_1
          
          
    cmd= "grep -A 28 'POSITION                                       TOTAL-FORCE' OUTCAR"
    #proc.wait()
    lines = call_subprocess(cmd)
    n_images = len(lines)/30
    lines = np.reshape(lines,(n_images,30))
    
    #print dir, "rs, miniAfter, total:",n_mini_rs, n_mini, n_images
    images=[]
    f = open('../../abnormal.out','a')
    traj = Trajectory('../../abnormal.traj','a')
    number = 0
    for items in lines:
        positions=[]
        forces=[]
        abnormal = False
        number +=1
        for i in range(0,30):
           if i >1 and i<17:
               fields = items[i].split()
               positions.append( [ float(fields[j]) for j in range(3) ] )
               forces.append( [ float(fields[j]) for j in range(3,6) ] )
           if i == 27:
              energy = float(items[27].split()[4])
              if energy < -100:
                 abnormal=True
                 f.write("%s %10d %20.6f\n"%(dir, number, energy))
        positions = np.array(positions)
        forces=np.array(forces)
        p = Atoms(elements, positions=positions)
        p.set_cell([[15.,0,0],[0,15.,0],[0,0,15.]],scale_atoms=False)
        p.set_pbc((True, True, True))
        calc = force_setter(energy=energy, forces=forces)
        p.set_calculator(calc)
        p.get_potential_energy()
        p.get_forces()
        images.append([p, energy, forces])
        if abnormal:
           traj.write(p)
   # path_images = images[n_mini_rs:n_images-n_mini]
    path_images=[]
    for config in dimer_configs:
       print dir,":", config
       for image in images[n_mini_rs:n_images-n_mini]:
          if rot_match(config, image[0], 0.01):
             #calc = force_setter(energy=energy, forces=forces)
             #config.set_calculator(calc)
             #config.get_potential_energy()
             #config.get_forces()
             path_images.append(image)
             break
    print len(path_images), len(dimer_configs)    
    rs_images = images[0:n_mini_rs]
    if n_mini_2 > 0:
       rs_images.extend(images[-n_mini:-n_mini_2])
       n_mini = n_mini_2
    if n_mini == 0:
       mini_images=[]
    else:
       mini_images=(images[-n_mini:])
    f.close()
    print dir, " rs, path, mini:", len(rs_images), len(path_images), len(mini_images) 
    return copy.deepcopy(rs_images), copy.deepcopy(path_images), copy.deepcopy(mini_images)

arg = sys.argv
current = os.getcwd()
cycle = int(arg[1])
jobid = int(arg[2])
current = os.getcwd()
#define elements
ref_atoms = read('pos.con', format='eon')
elements=ref_atoms.get_chemical_symbols()
state_main_dir = current+"/debug_results/"
os.chdir(state_main_dir)



job_listdir=[]
for f in glob.glob('*'):
   try:
      job_listdir.append(f)
   except:
      continue
job_listdir.sort()

rs_images=[]
path_images=[]
mini_images=[]
for dir in job_listdir:
   if int(dir.split('_')[0]) > cycle:
      break
   if int(dir.split('_')[1]) > jobid:
      continue
   os.chdir(state_main_dir+str(dir))
   rs, path, mini =read_atoms(dir, elements)
   rs_images.extend(rs)
   path_images.extend(path)
   mini_images.extend(mini)

from operator import itemgetter
rs_images    = sorted(rs_images, key=itemgetter(1))
path_images  = sorted(path_images, key=itemgetter(1))
mini_images  = sorted(mini_images, key=itemgetter(1))

os.chdir(current)
rs_traj = Trajectory('rs.traj','w')
path_traj = Trajectory('path.traj','w')
mini_traj = Trajectory('mini.traj','w')
for image in rs_images:
   rs_traj.write(image[0])
for image in path_images:
   path_traj.write(image[0])
for image in mini_images:
   mini_traj.write(image[0])
