#!/bin/sh

cwd="$(pwd)"
echo "PATH=\$PATH:$cwd/atoms_operator:$cwd/bin:$cwd/cp2k:$cwd/data_process:$cwd/jobscripts:$cwd/method_scripts:$cwd/vasp:$cwd/basinHopping_scripts:$cwd/crane:$cwd/expectra:$cwd/gnuplot_scripts:$cwd/neb_scripts" >>  ~/.bash_profile
echo "export PATH" >> ~/.bash_profile
