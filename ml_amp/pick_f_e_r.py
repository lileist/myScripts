from __future__ import division
import sys, random, copy
from ase.io import Trajectory
from ase import Atoms
import numpy as np

from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError


class forces_setter(Calculator):
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

def read_images(filename,state_number = None, mode = None):
    f = open(filename,'r')
    atoms = []
    ener = []
    fmins = []
    fmaxs = []
    fnorms = []
    cycle = -1
    while True:
        elements = []
        positions = []
        forces = []
        line = f.readline()
        if not line:
           break
        if cycle == -1:
           atom_numb = int(line.split()[0])
        line = f.readline()
        energy = float(line.split(":")[1])
        #read one particular structure assinged by state_number
        for i in range (atom_numb):
            line = f.readline()
            fields = line.split()
            elements.append(fields[0])
            positions.append( [ float(fields[j+1]) for j in range(3) ] )
            forces.append( [ float(fields[j+1]) for j in range(3,6) ] )
        elements = np.array(elements)
        positions = np.array(positions)
        ener.append(energy)
        p = Atoms(symbols=elements, positions=positions)
        calc = forces_setter(energy=energy, forces=np.array(forces))
        p.set_cell([[15.,0,0],[0,15.,0],[0,0,15.]],scale_atoms=False)
        p.set_pbc((True, True, True))
        p.set_calculator(calc)
        p.center()
        e  = p.get_potential_energy()
        fs = np.absolute(p.get_forces())
        fmin = np.amin(fs)
        fmax = np.amax(fs)
        fnorm = np.linalg.norm(np.sum(np.array(forces),axis=0))
        atoms.append([p, [energy, fnorm, fmin, fmax]])
        fmins.append(fmin)
        fmaxs.append(fmax)
        fnorms.append(fnorm)
        cycle += 1
    f.close()
    return atoms, np.array(ener), np.array(fnorms), np.array(fmins), np.array(fmaxs)

arg = sys.argv
configs, es, fnorms, fmins, fmaxs = read_images(arg[1])

print len(configs)
start_images = configs[np.argmin(es)]
start_coords = start_images[0].get_positions()
rs = []
#Tagged with distance to global minimum
for config in configs:
  print config[0]
  r = np.sum(np.linalg.norm(config[0].get_positions() - start_coords, axis = 1))
  config[1].extend([r])
  rs.append(r)

e_range = [np.amin(es), np.amax(es)]
fnorm_range = [np.amin(fnorms), np.amax(fnorms)]
fmin_range = [np.amin(fmins), np.amax(fmins)]
fmax_range = [np.amin(fmaxs), np.amax(fmaxs)]

#normalization
for config in configs:
   config[1] = [  config[1][0] ,
                ( config[1][1] - fnorm_range[0] ) / ( fnorm_range[1] - fnorm_range[0] ),
                ( config[1][2] - fmin_range[0] ) / ( fmin_range[1] - fmin_range[0] ),
                ( config[1][3] - fmax_range[0] ) / ( fmax_range[1] - fmax_range[0] ),
                  config[1][4]
               ]

duplicated_index = []
for i in range(len(configs)-1):
  config_1 = configs[i]
  print config_1[1]
  if i in duplicated_index:
     continue
  print i
  for j in range(i+1, len(configs)):
    config_2 = configs[j]
    if abs(config_2[1][0] - config_1[1][0]) < 0.02 and \
       abs(config_2[1][1] - config_1[1][1]) < 0.05 and \
       abs(config_2[1][2] - config_1[1][2]) < 0.05 and \
       abs(config_2[1][3] - config_1[1][3]) < 0.05 and \
       abs(config_2[1][4] - config_1[1][4]) < 0.05:
       duplicated_index.append(j)

print len(configs)
for index in sorted(duplicated_index, reverse = True):
   del configs[index]
print len(configs)
traj = Trajectory('MD_images.traj', 'w')
#test_traj = Trajectory('test_images.traj', 'w')
#del configs[:min_index]

"""
random.shuffle(configs_list)
n_split = int(0.8 * (len(configs_list) + len(high_f))) - len(high_f)
images_train = configs_list[:n_split]
images_test = configs_list[n_split:]
images_train.extend(high_f)
"""
for image in configs:
    #i+=1
    #train.write("%10d %12.8f\n"%(i, image.get_potential_energy()))
    traj.write(image[0])
