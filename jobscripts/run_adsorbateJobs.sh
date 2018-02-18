#!/bin/sh
set -e

#if [ -e runs0 ]; then
#    rm -rf runs0
#fi
#mkdir -p runs0
#cd runs0
for i in {1..6}
do
    file=$(printf "POSCAR%i" $i)
    mkdir -p run-$i
    /bin/echo -n run-$i
    cp $file run-$i/POSCAR
    cp INCAR run-$i
    cp POTCAR run-$i
    cp KPOINTS run-$i
    cp qsub.hf run-$i
    cd run-$i
    qsub qsub.hf
    cd ..
done
