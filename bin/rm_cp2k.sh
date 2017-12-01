#!/bin/sh
find ./ -name "*Hessian"|xargs rm
find ./ -name "*bak*"|xargs rm
find ./ -name "*WFN*"|xargs rm
find ./ -name "*wfn*"|xargs rm
find ./ -name "machines*"|xargs rm
#find ./ -name "jobid"|xargs rm
find ./ -name "*restart"|xargs rm
find ./ -name "*cube*"|xargs rm
#find ./ -name "suppl.out"|xargs rm
find ./ -name "*-vel-1.xyz"|xargs rm
