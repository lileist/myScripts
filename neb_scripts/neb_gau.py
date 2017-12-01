#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator.

http://www.cp2k.org
Author: Ole Schuett <ole.schuett@mat.ethz.ch>
"""

from __future__ import division, print_function
import os
import sys
from ase.test import NotAvailable
from ase import atoms
from ase import units
from ase.structure import molecule
from ase.calculators.gaussian import Gaussian
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from tsase import neb
"""
neb_gau [reactantFile] [productFile] [#ofCores] [restart]
"""

def main():
    arg = sys.argv
    # Basically, the entire CP2K input is passed in explicitly.
    # Disable ASE's input generation by setting everything to None.
    # ASE should only add the CELL and the COORD section.
    calc=Gaussian(method='M062x',
                basis='6-311++G(d,p)',
                nproc = int(arg[3]),
                charge = 0,
                force = 'force',
                chk = '3w.chk',
                label='3w')
    p1 = read(arg[1],index=0,format='xyz')
    #p1.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False,fix=None)
    #p1.set_pbc((True, True, True))
    p1.set_calculator(calc)
    p2 = read(arg[2],index=0,format='xyz')
    #p2.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False,fix=None)
    #p2.set_pbc((True, True, True))
    p2.set_calculator(calc)
#    MaxwellBoltzmannDistribution(slab_from_file, 0.5 * 300 * units.kB, force_temp=True)
#    energy_start = atoms.get_potential_energy() + atoms.get_kinetic_energy()
#    dyn = VelocityVerlet(slab_from_file, 0.5 * units.fs)
    #def print_md():
    #    energy = atoms.get_potential_energy() + atoms.get_kinetic_energy()
    #    print("MD total-energy: %.10feV" %  energy)
    #dyn.attach(print_md, interval=1)
#    dyn.run(20)
    nim = 9  # number of images, including end points
    band = neb.ssneb(p1, p2, numImages = nim, method = 'normal', ss=False)
    if len(arg)>4:
       restart = True
    else:
       restart = False
# to restart, uncomment the following lines which read the previous optimized images into the band
    if restart:
       for i in range(1,nim-1):
           filename = str(i)+'.CON'
           b = read(filename,format='vasp')
           band.path[i].set_positions(b.get_positions())
           band.path[i].set_cell(b.get_cell())

#    opt = neb.qm_ssneb(band, maxmove =0.2, dt=0.05)
    opt = neb.fire_ssneb(band, maxmove =0.2, dtmax = 0.1, dt=0.1)
    opt.minimize(forceConverged=0.1, maxIterations = 1000)

main()
# EOF

