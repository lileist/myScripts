#!/bin/sh

for i in 1 2 3 4 5
do
cd 0$i
cp CONTCAR POSCAR
cd ..
echo -e $i
done
