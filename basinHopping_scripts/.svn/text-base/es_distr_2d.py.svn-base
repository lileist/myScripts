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
    energy = []
    s =[]
    for line in f:
        if line.startswith('#'):
           continue
        fields = [ field for field in line.split()]
        energy.append(float(fields[6]))
        s.append(float(fields[7]))
        #dot = [float(fields[6]), float(fields[7])]
    f.close()
    return energy, s

def printMatrix(a, e_min, e_max, s_min, s_max):
   print "# Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]"
   print "# "+ ("%15.6f" % e_min) +("  %15.6f" % e_max) +("  %15.6f" % s_min) +("  %15.6f" % s_max)
   rows = a.shape[0]
   cols = a.shape[1]
   for i in range(0,rows):
      for j in range(0,cols):
         print("%8.5f" %a[i,j]),
      print
   print  

def main():
    arg = sys.argv
    inputfile = arg[1]
    de = float(arg[2])
    ds = float(arg[3])
    energy,s = read_dots(inputfile)
    e_min = min(energy)
    e_max = max(energy)
    s_min = min(s)
    s_max = max(s)

    #initialize e_gr and s_gr
    w = int((e_max - e_min)/de)+1
    h = int((s_max - s_min)/ds)+1
    es_gr = [[0 for x in range(w)] for y in range (h)]
    #analyze
    for i in xrange (len(energy)):
        for k in xrange (int((s_max - s_min)/ds)+1):
            for j in xrange (int((e_max - e_min)/de)+1):
                if energy[i] >= e_min + float(j)*de and energy[i] < e_min + float((j+1))*de \
                   and s[i] >= s_min + float(k)*ds and s[i] < s_min+float((k+1))*ds:
                   es_gr[k][j] += 1
    #normalize
    for i in xrange (int((e_max - e_min)/de)+1):
        for j in xrange (int((s_max - s_min)/ds)+1):
            es_gr[j][i] = float(es_gr[j][i])/float(len(energy))
    a=np.array(es_gr)
    printMatrix(a, e_min, e_max, s_min, s_max)

if __name__ == '__main__':
    main()

        
