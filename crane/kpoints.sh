#!/bin/sh

for i in 6 8 10 12 14 16 18 20 22 24 26
do
rm WAVECAR
#i=$(echo "$i"|bc -l) 
#a=$(echo "$i"|bc -l) 
## echo "$a $b"


cat >KPOINTS << CATY
K-Points
0
Gamma
$i $i 1
0 0 0
CATY
echo '==========poscar============'  >> out
cat POSCAR  >> out
echo '============================='  >> out

mpirun -np 4 ~/bin/vasp5.3.3 </dev/null >> out

echo -n "0 0 $i :  " >> volum_energy.txt
grep "energy without entropy" OUTCAR |tail -1|awk '{print $8}' >> volum_energy.txt
#grep "energy without entropy" OUTCAR |tail -1
echo $i >> contcar-save
cat CONTCAR > $i
echo ' '
grep "energy without entropy" OUTCAR |tail -1   >> out
grep "EENTRO =" OUTCAR | tail -1   >> out
done

