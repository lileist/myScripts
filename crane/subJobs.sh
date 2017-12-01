#!/bin/sh

for i in 2 3 4 5 6 7
do
cd $i
sbatch run-vasp-intel.slurm
cd ..
echo -e $i
done
