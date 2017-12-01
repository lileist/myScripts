#!/usr/bin/env python
import sys
sys.path.append("../../")
import os

#from ase.test import NotAvailable
from ase.io import read
from ase.units import kB
from ase.optimize.lbfgs import LBFGS
from ase.calculators.emt import EMT
from ase.atoms import Atoms
from expectra.cal_exafs import Expectra
#from expectra.io import read_xdatcar, read_con, read_chi
from ase.utils.geometry import sort
#from expectra.paretoOPT import ParetoLineOptimize
from expectra.basin_surface import BasinHopping
from tsase.calculators.lammps_ext import LAMMPS

def log_atoms(atoms, log_trajectory):
    log_trajectory.write("%d\n" % (len(atoms)))
    log_trajectory.write(" \n")
    for atom in atoms:
        log_trajectory.write("%s  %15.6f  %15.6f  %15.6f\n" % (atom.symbol,
                                 atom.x, atom.y, atom.z))
    log_trajectory.flush()
def main():
    #read in geometries
    p1 = read('CONTCAR',index=0,format='vasp')

#    charges = [setcharge(a.symbol) for a in p1]
#    p1.set_initial_charges(charges)
#    potfile=['library.meam','Au-Rh.meam']
    #calculator for geometry optimization and MD
     
    bh = BasinHopping(atoms=p1,
                      ncore=2,
                      specorder = ['Au', 'Rh'],
                      switch = True,
                      active_ratio = 0.08,
                      #jumpmax = 21,
                      #adjust_step_size = 10,
                      #significant_structure = True,
                      elements_lib = ['Au', 'Rh'])
    #run job
    log_trajectory = open('test.xyz', 'w')
    for i in range (10):
        p1.set_chemical_symbols(bh.random_swap(p1.get_chemical_symbols()))
        log_atoms(sort(p1), log_trajectory)
if __name__ == '__main__':
    main()
    
