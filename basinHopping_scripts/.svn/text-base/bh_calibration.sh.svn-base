#!/bin/sh
set -e

cd runs0
for i in {0..4}
do
    file=$(printf "%i.xyz" $i)
    mkdir -p run-$i
    /bin/echo -n run-$i
    cp ../../../random_stru/$file run-$i/pos.xyz
    cp ../EMT_expectra.py run-$i
    cp ../basin_surface.py run-$i
    cp ../PdAu_sc.eam.fs run-$i
    cd run-$i
    ./EMT_expectra.py > stdout.dat
    #serialsub eonclient
    #tail -n 1 pot_log.dat >> ../../forcecalls.out
    cd ..
done
