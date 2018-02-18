#!/usr/bin/env python
"""
This code is used to analyze the distribution of energy and s (chi_deviation).
How to use:
    gr_ES [inputfile_name] [de] [ds]
"""
import sys
import os
import numpy as np

def read_dots(start,end,increment):
    energy = []
    s =[]
    directory = os.getcwd()
    for i in range(start,end, increment):
       bh_file=directory+"/run-"+str(i)+'/pot'
       if os.path.exists(bh_file):
          f_bh = open(bh_file+'/pot_log.dat','r')
       else:
          #num_bh[str(i)] = 0
          print 'run-',i, 'is not run'
          continue
       lines = f_bh.readlines()
       for line in lines:
           if line.startswith('#'):
              continue
           fields = [ field for field in line.split()]
           bh_step = fields[1]
           energy.append(float(fields[6]))
           s.append(float(fields[7]))
           if int(bh_step) == 20000:
              break
       f_bh.close()
    return np.array(energy), np.array(s)

def printMatrix(a, e_min, e_max, e_mean, s_min, s_max,s_mean):
   print "# Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]"
   print "# "+ ("%15.6f" % e_min) +("  %15.6f" % e_max) +("  %15.6f" % s_min) +("  %15.6f" % s_max) +("  %15.6f" % e_mean)+("  %15.6f" % s_mean)
   rows = a.shape[0]
   cols = a.shape[1]
   for i in range(0,rows):
      for j in range(0,cols):
         print("%8.5f" %a[i,j]),
      print
   print  

def main():
    arg = sys.argv
    de = float(arg[3])
    ds = float(arg[4])
    increment = 5
    energy,s = read_dots(int(arg[1]),int(arg[2]),increment)
    e_min = min(energy)
    e_max = max(energy)
    s_min = min(s)
    s_max = max(s)
    e_mean = np.mean(energy)
    s_mean = np.mean(s)
    #initialize e_gr and s_gr
    w = int((e_max - e_min)/de)+1
    h = int((s_max - s_min)/ds)+1
    es_gr = [[0 for x in range(w)] for y in range (h)]
    #analyze
    t_numb=0
    for i in xrange (len(energy)):
        for k in xrange (int((s_max - s_min)/ds)+1):
            for j in xrange (int((e_max - e_min)/de)+1):
                if energy[i] >= e_min + float(j)*de and energy[i] < e_min + float((j+1))*de \
                   and s[i] >= s_min + float(k)*ds and s[i] < s_min+float((k+1))*ds:
                   es_gr[k][j] += 1
                   t_numb+=1
    print t_numb
    #normalize
    for i in xrange (int((e_max - e_min)/de)+1):
        for j in xrange (int((s_max - s_min)/ds)+1):
            es_gr[j][i] = float(es_gr[j][i])/float(len(energy))
    a=np.array(es_gr)
    printMatrix(a, e_min, e_max, e_mean, s_min, s_max, s_mean)

if __name__ == '__main__':
    main()

        
