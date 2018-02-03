#qstat -u "*" |grep 'au147_fs'|awk -F"\t" '{if ($1) print $1}' >out
import sys, os
arg = sys.argv

f = open(arg[1],'r')

lines = f.readlines()
for line in lines:
   print line
   os.system("qdel "+line.split()[0])

