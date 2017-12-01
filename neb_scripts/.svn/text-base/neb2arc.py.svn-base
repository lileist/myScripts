#!/usr/bin/env python
import sys
import os
import shutil
from ase.io import Trajectory
from ase.io import read
from ase import Atoms

""" convert NEB results to a trajectory file (.arc) which can be read by MS
    Note: arc file has very strict format requirement
          i.e., xyz coordinates cannot be too forward or backward;
                'XXXX 1' have to start from column 52
"""

def log_atoms(atoms):
    outfile = open('movie.arc', 'w')
    for i in range(len(atoms)):
        if i == 0:
           outfile.write("!BIOSYM archive 3\n")
           outfile.write("PBC=ON\n")
        outfile.write(" \n")
        outfile.write("!DATE\n")
        #TODO: output non-orthogonal cells
        lattice = atoms[i].cell
        outfile.write("PBC  %8.4f  %8.4f  %8.4f  %9.4f  %9.4f  %9.4f\n" %
                      (lattice[0][0], 
                       lattice[1][1],
                       lattice[2][2],
                       90.000, 90.000, 90.000))
        for atom in atoms[i]:
            x = atom.x
            y = atom.y
            z = atom.z
            if abs(atom.x - atoms[0][atom.index].x) > lattice[0][0]/2:
               if atom.x > atoms[0][atom.index].x:
                  x = atom.x - lattice[0][0]
               else:
                  x = atom.x + lattice[0][0]
            if abs(atom.y - atoms[0][atom.index].y) > lattice[1][1]/2:
               if atom.y > atoms[0][atom.index].y:
                  y = atom.y - lattice[1][1]
               else:
                  y = atom.y + lattice[1][1]
            if abs(atom.z - atoms[0][atom.index].z) > lattice[2][2]/2:
               if atom.z > atoms[0][atom.index].z:
                  z = atom.z - lattice[2][2]
               else:
                  z = atom.z + lattice[2][2]
            x = round(x,8)
            y = round(y,8)
            z = round(z,8)
            coordinates = atom.symbol + '      '+str(x).ljust(12)+'    '+str(y).ljust(12)+'    '+str(z).ljust(12)
            coordinates = coordinates.ljust(51) + 'XXXX 1      xx      ' + atom.symbol +'  '+ '0.000'
            outfile.write("%s\n" % (coordinates))
            #outfile.write("%s      %12.8f    %12.8f    %12.8f XXXX 1      xx      %s  0.000\n" % 
            #               (atom.symbol, atom.x, atom.y, atom.z, atom.symbol))
        outfile.write("end\n")
        outfile.write("end\n")
        
def main():
    arg = sys.argv
    num_image = int(arg[1])
    cwd = os.getcwd()
    if num_image is None:
        print "num_image is empty"
    atoms = []
    if len(arg)>3:
       vasp = False
    else:
       vasp = True
    if vasp:
       for image in range(0, num_image):
            if image==0 or image == num_image-1:
               target_file = 'POSCAR'
            else:
               target_file = arg[2]
            if image > 9:
               cwd_1 = cwd + '/'+str(image)+'/'
            else:
               cwd_1 = cwd + '/0'+str(image)+'/'
            src = cwd_1+target_file
            dst = cwd +'/'+arg[2]
            shutil.copy(src, dst)
            print(src+' --> '+dst)
            p1 = read(filename=arg[2], index=0, format = 'vasp')
            atoms.append(p1)
    else:
       for image in range(0, num_image+2):
           if image==0:
              p1 = read(filename=arg[2], index=0, format = 'xyz')
           elif image==num_image+1:
              p1 = read(filename=arg[3], index=0, format = 'xyz')
           else:
              p1 = read(filename=str(image)+'.CON', index=0, format = 'vasp')
           atoms.append(p1)
    log_atoms(atoms)
if __name__ == '__main__':
    main()
