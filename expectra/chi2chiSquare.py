#!/usr/bin/env python
"""
This code is used to extract geometry and exafs for given energy and s (chi_deviation).
How to use:
   dot_selection  [energy] [s] [pl_steps] [bh_steps] [species_number]
"""
import sys
import os
import numpy
from expectra.io import read_chi
from expectra.cal_exafs import save_result

def main():
    arg = sys.argv
    inputfile = arg[1]
    k, chi = read_chi(inputfile)
    reverse = int(arg[2])
    if reverse == 1:
       save_result(k, numpy.divide(chi, numpy.square(k)), arg[1]+'_1')
    else:
       save_result(k, numpy.multiply(numpy.square(k), chi), arg[1]+'2')

if __name__ == '__main__':
    main()

        
