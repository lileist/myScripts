#!/bin/sh
set -e

#if [ -e runs0 ]; then
#    rm -rf runs0
#fi
#mkdir -p runs0
#cd runs0
for i in {0..9}
do
    file=$(printf "%i.xyz" $i)
    mkdir -p run-$i
    /bin/echo -n run-$i
    cp ../random_stru/$file run-$i/pos.xyz
    cp bh.py run-$i
    cp qsub.hf run-$i
    cp basin.py run-$i
    cp PdAu_sc.eam.fs run-$i
    cp exp_exafs.chi run-$i
    cd run-$i
    qsub qsub.hf
    #serialsub eonclient
    #tail -n 1 pot_log.dat >> ../../forcecalls.out
    cd ..
done
