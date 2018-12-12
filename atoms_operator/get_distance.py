
from ase.io import read
import sys
import numpy
args = sys.argv
atoms = read(args[1],index=":")

distance = []

for image in atoms:
   distance.append(image.get_distance(int(args[2]), int(args[3])))
print "#averge:",numpy.mean(numpy.array(distance))
for i in range(len(distance)):
   print i, distance[i]

