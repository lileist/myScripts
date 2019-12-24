import sys
from expectra.atoms_operator import match, single_atom
from ase.io import Trajectory, read
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError
import numpy as np
import argparse

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
        p = Atoms(symbols=elements, positions=positions)
        calc = forces_setter(energy=energy, forces=np.array(forces))
        p.set_cell([[20.,0,0],[0,20.,0],[0,0,20.]],scale_atoms=False)
        p.set_pbc((True, True, True))
        p.set_calculator(calc)
        p.center()
        e  = p.get_potential_energy()
        fs = np.absolute(p.get_forces())
        fmin = np.amin(fs)
        fmax = np.amax(fs)
        if fmax > 10.0:
           continue
        ener.append(energy)
        fnorm = np.linalg.norm(np.sum(np.array(forces),axis=0))
        atoms.append([p, [energy, fnorm, fmin, fmax]])
        fmins.append(fmin)
        fmaxs.append(fmax)
        fnorms.append(fnorm)
        cycle += 1
    f.close()
    return atoms, np.array(ener), np.array(fnorms), np.array(fmins), np.array(fmaxs)

def read_traj(filename):
    trajs = read(filename, index=":")
    #trajs = Trajectory(filename, "r")
    ener = []
    atoms=[]
    fnorms = []
    fmins = []
    fmaxs = []
    for p in trajs:
        #p.set_cell([[20.,0,0],[0,20.,0],[0,0,20.]],scale_atoms=False)
        #p.set_pbc((True, True, True))
        #p.center()
        e  = p.get_potential_energy()
        forces = p.get_forces()
        fs = np.absolute(forces)
        fmin = np.amin(fs)
        fmax = np.amax(fs)
        if fmax > 10.0:
           continue
        ener.append(e)
        fnorm = np.linalg.norm(np.sum(forces,axis=0))
        atoms.append([p, e])
        fmins.append(fmin)
        fmaxs.append(fmax)
        fnorms.append(fnorm)
    return atoms, np.array(ener), np.array(fnorms), np.array(fmins), np.array(fmaxs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--traj', type=str, nargs='+', metavar='trajFile', 
                        default=None,
                        help='filename of the trajectory')
    parser.add_argument('--runtype', type=str, nargs='+', metavar='type', 
                        default=None,
                        help='run type: split or match')
    args = parser.parse_args()
    print(args)
    if args.traj[0].split('.')[1]=='xyz':
       configs, es, fnorms, fmins, fmaxs = read_images(args.traj[0])
    if args.traj[0].split('.')[1]=='traj':
       configs, es, fnorms, fmins, fmaxs = read_traj(args.traj[0])
    if args.runtype[0] == 'split':
       configs.sort(key=lambda x:x[1])
       dn = int(len(configs)/14)
       for i in range(14):
          print(i)
          test_traj = Trajectory(str(i)+'.traj','w')
          start = dn*i
          if i == 6:
             end=len(configs)
          else:
             end   = dn*(i+1)
          for config in configs[start : end]:
             test_traj.write(config[0])

    if args.runtype[0] == 'match':
       outtraj = Trajectory('filter.traj','w')
       for i in range(len(configs)-1):  
       #calc = force_setter(energy=config.get_potential_energy(), forces=config.get_forces())
       #config.set_cell([[40.,0,0],[0,40.,0],[0,0,40.]],scale_atoms=False)
       #config.set_pbc((True, True, True))
       #config.center()
       #config.set_calculator(calc)
          e0 =configs[i][0].get_potential_energy()
          fs0=configs[i][0].get_forces()
          matched = False
          for j in range(i+1, len(configs)):
             e1 =configs[j][0].get_potential_energy()
             if abs(e1-e0) <= 0.05:
                matched=match(configs[i][0], configs[j][0], 0.1, 3.3, indistinguishable=True)
                if matched:
                   break
          print(i, matched)
          if matched:
             continue
          outtraj.write(configs[i][0])
    #if matched:
    #   print "matched"
    #print matched
if __name__ == '__main__':
    main()


