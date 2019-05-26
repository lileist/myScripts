#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Test suit for the CP2K ASE calulator with neb calculator from TSASE

http://www.cp2k.org
Author: Ole Schuett <ole.schuett@mat.ethz.ch>
"""

from __future__ import division, print_function
import os

from ase.test import NotAvailable
from ase import atoms
from ase import units
from ase.structure import molecule
from amp import Amp
from ase.io import read
from tsase import neb


def main():
    calc = Amp.load('amp.amp')
    """Read the reactant coordinates"""
    p1 = read('POSCAR_rs',index=0,format='vasp')
#    p1.set_cell([[20,0,0],[0,20,0],[0,0,20]],scale_atoms=False)
#    p1.set_pbc((True, True, True))
    p1.set_calculator(calc)

    """Read the product coordinates"""
    p2 = read('POSCAR_fs',index=0,format='vasp')
    p2.set_calculator(calc)

    nim = 12  # number of images, including end points
    band = neb.ssneb(p1, p2, numImages = nim, method = 'normal', ss=False)

# to restart, uncomment the following lines which read the previous optimized images into the band
#    for i in range(1,nim-1):
#        filename = str(i)+'.CON'
#        b = read(filename,format='vasp')
#        band.path[i].set_positions(b.get_positions())
#        band.path[i].set_cell(b.get_cell())

#    opt = neb.qm_ssneb(band, maxmove =0.1, dt=0.05)
#Uncomment to use fire optimization algorithm
    opt = neb.fire_ssneb(band, maxmove =0.1, dtmax = 0.05, dt=0.05)
    opt.minimize(forceConverged=0.10, maxIterations = 400)


main()
# EOF

