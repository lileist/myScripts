#!/bin/sh
for i in {0..39}
do
   cd run-$i
   #rm au55* new* out pareto* *lammps*
   mv pot/pot* .
   rm -rf basin_surface.py geo* basin_surface.pyc  configs/  data.lmp.initialize  EMT_expectra.py*  exafs/  in.lmp.initialize  lowest.traj  PdAu_sc.eam.fs  pos.xyz time_log.dat pot err  ll_out  out
   cd ..
done
