import numpy as np
import os, sys

def main():
    directory = os.getcwd()
    arg = sys.argv
    output = open("pots.dat",'w')
    output.write("step    energy       chi        pseudoPot    Umin\n")
    increment = 5
    for i in range(int(arg[1]), int(arg[2]), increment):
       bh_file=directory+"/run-"+str(i)+'/pot'
       print bh_file
       if os.path.exists(bh_file):
          pots = np.loadtxt(bh_file+'/pot_log.dat',skiprows=1,usecols=(6,7,9,10,16))
          #pots = np.genfromtxt(bh_file+'/pot_log.dat',delimiter=' ',skip_header=1,usecols=(6,7,9,10,16))
       else:
          print 'run-',i, 'is not run'
          continue
       for pot in pots:
          output.write("%15.6f %15.6f %15.6f %15.6f %15.6f\n"%(pot[0],pot[1],pot[2],pot[3],pot[4]))
if __name__ == '__main__':
    main()
    
