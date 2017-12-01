#!/bin/sh

for i in 1 2 3 4 5
do
cd 0$i
grep FORCES: OUTCAR >forces
tail -1000 OUTCAR >outcar
rm OUTCAR
cd ..
echo -e $i
done
