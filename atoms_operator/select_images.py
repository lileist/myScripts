
from ase import Atoms
import sys
from set_forces import force_setter

def read_images(filename,state_number = None, mode = None):
    f = open(filename,'r')
    atoms = []
    ener = []
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
        ener.append(float(line.split()[0]))
        #read one particular structure assinged by state_number
        if state_number is not None:
           if cycle == state_number:
              for i in range (atom_numb):
                  line = f.readline()
                  fields = line.split()
                  elements.append(fields[0])
                  positions.append( [ float(fields[j+1]) for j in range(3) ] )
              elements = np.array(elements)
              positions = np.array(positions)
              return Atoms(elements, positions=positions)
           else:
              for i in range (atom_numb):
                  f.readline()
        else:
           for i in range (atom_numb):
               line = f.readline()
               fields = line.split()
               elements.append(fields[0])
               positions.append( [ float(fields[j+1]) for j in range(3) ] )
               forces.append( [ float(fields[j+1]) for j in range(3,6) ] )
           elements = np.array(elements)
           positions = np.array(positions)
           p = Atoms(elements, positions=positions)
           calc = force_setter()
           p.set_calculator(calc)
           print p.get_potential_energy()
           print p.get_forces()
           atoms.append(p)
        cycle += 1
    f.close()
    return atoms

