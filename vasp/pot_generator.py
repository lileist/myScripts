#!/usr/bin/env python
import sys
def main():
  arg = sys.argv
  dir = arg[1]
  species = arg[2]
#  spin = int(arg[2])
  pot_input = open('LOCPOT', 'r')
  while True:
      line = pot_input.readline()
      if not line:
         break
      for i in range(6):
          line = pot_input.readline()
      ions = [int(field) for field in line.split()]
      print "ions:", ions
      n_ions = sum(ions)
      print "total number of ions:", n_ions
      line = pot_input.readline()
      for i in range(n_ions):
          pot_input.readline()
      pot_input.readline()
      vlocal = []
      vlocal.append(0.0)
#      for i in range(spin):
      fields = [int(field) for field in pot_input.readline().split()]
      ngx = fields[0]
      ngy = fields[1]
      ngz = fields[2]
      print "ngx ngy ngz:",fields 
      nplwv = ngx*ngy*ngz
      #first element is not used
      for j in xrange(int(nplwv/5)+1):
          vlocal.extend([float(field) for field in pot_input.readline().split()])
      print len(vlocal)
      if len(vlocal)>=nplwv+1:
         break
  if idir == 'x':
     nout = ngx
  elif idir == 'y':
     nout = ngy
  else:
     nout = ngz
  scale = 1./float(nplwv/nout)
  vav = []
  vav.append(0.0)
  for i in xrange(nout):
     vav.append(0.0)
  if idir == 'x':
    for i in xrange(1, ngx+1):
        for j in xrange(1, ngy+1):
            for k in xrange(1, ngz+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[i] = vav[i]+ vlocal[ipl] * scale
  elif idir == 'y':
    for j in xrange(1, ngy+1):
        for k in xrange(1, ngz+1):
            for i in xrange(1, ngx+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[j] = vav[j]+ vlocal[ipl] * scale
  elif idir == 'z':
    for k in xrange(1, ngz+1):
        for j in xrange(1, ngy+1):
            for i in xrange(1, ngx+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[k] = vav[k]+ vlocal[ipl] * scale
  output = open('loc_pot.dat','w')
  for i in xrange(1, len(vav)):
      output.write("%d %15.6f\n"%(i, vav[i]))
if __name__ == '__main__':
    main()
