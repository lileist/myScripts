import numpy as np
import sys
arg = sys.argv
#data = np.genfromtxt(arg[1],delimiter=None,skip_header=2)
f = open(arg[1],'r')
data = []
while True:
   line = f.readline()
   if not line:
      break
   if "min:" not in line:
     continue
   data.append(float(line.split()[10]))
print np.mean(data[300:501])
