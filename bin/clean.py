#!/usr/bin/env python
"""
This code is used to analyze the distribution of energy and s (chi_deviation).
How to use:
    gr_ES [inputfile_name] [de] [ds]
"""
import sys
import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def rm_exafs(mode=None):
    os.system('rm ll_out opt_lammps* out md_lammps* new_au55_all.xyz.gz new_au55_exafs_all.dat.gz time_log* paretoAtoms.traj paretoLine.dat blz_sample.dat')
    os.system('rm -rf configs exafs visited_configs.dat trj_lammps new_configs* data_lammps log.lammps chi.dat au55_all.xyz au55_exafs_all.dat')
    os.system('rm geo_opt.log.gz')
def main():
    arg = sys.argv
    filetype = arg[1]
    for i in range(int(arg[2]), int(arg[3])):
       with cd('run-'+str(i)):
            if filetype == 'exafs':
               rm_exafs()
if __name__ == '__main__':
    main()

        
