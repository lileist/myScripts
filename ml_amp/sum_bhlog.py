
bh_log = []
for i in range(14):
  f = open('out_'+str(i), 'r')
  lines = f.readlines()
  print(i,':', len(lines))
  if i ==0:
     for j in range(len(lines)):
       bh_log.append(float(lines[j].rstrip().split(' ')[-1]))
  else:
     for j in range(len(lines)):
       bh_log[j]+=float(lines[j].rstrip().split(' ')[-1])
  f.close()
out = open('bh_log.dat', 'w')
for i in range(len(bh_log)):
  out.write("{:6d} {:12.6f}\n".format(i, bh_log[i]))
