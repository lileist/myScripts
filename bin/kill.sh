#!/bin/sh
set -e

#if [ -e runs0 ]; then
#    rm -rf runs0
#fi
#mkdir -p runs0
#cd runs0
for i in {1348028..1348067}
do
    qdel $i
    #serialsub eonclient
    #tail -n 1 pot_log.dat >> ../../forcecalls.out
done
