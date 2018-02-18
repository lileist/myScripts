import sys
import os
import subprocess

arg = sys.argv
flag = 'reached'
for i in range(int(arg[1])):
   cmd = 'run-'+str(i)
   grep_flag = 'grep '+flag+' '+arg[2]
   os.system(cmd)
   proc = subprocess.Popen(grep_flag,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
   output = proc.communicate()
   #lines = output[0].split('\n')
   if 'reached required accuracy' in output[0]:
      print 'Job '+cmd, 'is done'
      os.system('cd ..')
      continue
   else:
      os.system('cp CONTCAR POSCAR')
      os.system('qsub qsub.hf')
      print 'Rerun ', cmd
      os.system('cd ..')

