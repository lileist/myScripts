
import numpy
import sys, os

arg = sys.argv
f = open('u_log.dat', 'r')
lines = f.readlines()
nelect = []
pot = []
energy = []
for line in lines:
    if line.startswith('step'):
       continue
    fields = line.split()
    nelect.append(float(fields[1]))
    pot.append(float(fields[2]))
    energy.append(float(fields[5]))
fit_paras = numpy.polyfit(nelect, pot, 1)
k = fit_paras[0]
b = fit_paras[1]
print (float(arg[1]) - b)/k
fit_paras_1 = numpy.polyfit(pot, energy, 1)
k = fit_paras_1[0]
b = fit_paras_1[1]
print k*float(arg[1])+b
