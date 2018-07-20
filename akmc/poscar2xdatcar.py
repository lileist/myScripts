#used to convert .poscar (generated from eon -m) to XDATCAR
# inputtrajectory  n_atom  n_traj
import sys

arg = sys.argv
f = open(arg[1],'r')

fo = open('XDATCAR','w')

n_atom = int(arg[2])
for n_traj in range(int(arg[3])):
  for i in range(8):
     line = f.readline()
     if n_traj >0:
        continue
     if i == 2:
       l = float(line.split()[0])
     if i ==6:
       fo.write("%s %d\n"%("Direct configuration=", n_traj+1))
       continue
     if i >6:
       continue
     fo.write(line) 
  if n_traj >0:
    fo.write("%s %d\n"%("Direct configuration=", n_traj+1))
  for i in range(n_atom):
     fields = f.readline().split()
     x = float(fields[0])/l
     y = float(fields[1])/l
     z = float(fields[2])/l
     fo.write("%12.8f   %12.8f   %12.8f\n"%(x,y,z))
  
  
