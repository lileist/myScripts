import os,sys

args = sys.argv
f = open('config.ini','r')
lines = f.readlines()
for i in range(int(args[1])):
  os.mkdir('opt-'+str(i))
  fout = open('opt-'+str(i)+'/config.ini','w')
  for line in lines:
     if 'structure_slice' in line:
        fout.write("{:s}\n".format('structure_slice = '+str(i*100)+':'+str((i+1)*100)+':1'))
        continue
     fout.write(line)
