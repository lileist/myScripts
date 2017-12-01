#!/bin/sh
set -e

#if [ -e runs0 ]; then
#    rm -rf runs0
#fi
#mkdir -p runs0
#cd runs0
for i in {0..17}
do
    file=$(printf "%i.xyz" $i)
    rm run-$i/geo_opt
done
