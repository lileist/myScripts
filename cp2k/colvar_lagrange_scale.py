#!/usr/bin/env python
"""
This code is used to analyze the distribution of energy and s (chi_deviation).
How to use:
    gr_ES [inputfile_name] [de] [ds]
"""
import sys
import os
import numpy as np

def read_dots(filename):
    f = open(filename)
    cv_auxi = []
    cv =[]
    for line in f:
        if line.startswith('#'):
           continue
        fields = [ field for field in line.split()]
        cv_auxi.append(float(fields[1]))
        cv.append(float(fields[2]))
    f.close()
    return np.array(cv_auxi), np.array(cv)


def main():
    arg = sys.argv
    inputfile = arg[1]
    temp = float(arg[2])/315775.0
    k = float(arg[3])
    cv_auxi,cv = read_dots(inputfile)
    
    cv_vari = np.average(np.square(cv_auxi-cv))
    cv_target = (np.average(np.square(cv_auxi))-(np.average(cv_auxi))**2)*temp/k
    print cv_vari, cv_target, cv_target*k/cv_vari
if __name__ == '__main__':
    main()

        
