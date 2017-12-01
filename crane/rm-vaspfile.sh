#!/bin/sh
find ./ -name "DOS*"|xargs rm
find ./ -name "PRO*"|xargs rm
find ./ -name "vasprun*"|xargs rm
find ./ -name "CHG*"|xargs rm
find ./ -name "WAVECAR"|xargs rm
find ./ -name "XDATCAR"|xargs rm
find ./ -name "OUTCAR"|xargs rm
find ./ -name "EIGEN*"|xargs rm
find ./ -name "IBZ*"|xargs rm
find ./ -name "OSZ*"|xargs rm
find ./ -name "PCD*"|xargs rm
find ./ -name "out"|xargs rm
find ./ -name "machine*"|xargs rm
find ./ -name "jobid"|xargs rm
