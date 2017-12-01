#!/usr/bin/env python
"""
This code is used to extract geometry and exafs for given energy and s (chi_deviation).
How to use:
   dot_selection  [energy] [s] [pl_steps] [bh_steps] [species_number]
"""
import sys
import os

def read_dots(filename, energy, s):
    f = open(filename)
    #print energy," ", s
    for line in f:
        if line.startswith('#'):
           continue
        fields = [ field for field in line.split()]
        #print fields[6], " ", fields[7]
        if fields[6] == energy and fields[7] == s:
           return fields[1]
    f.close()
    return None

def read_atoms(filename, step, node, bh_step):
    f = open(filename, 'r')
    atom_numb = None
    while True:
        line = f.readline()
        if not line:
           break
        fields = [ field for field in line.split()]
        if len(fields)==1 and atom_numb is None:
           atom_numb = int(fields[0])

        if 'node: None  step:' in line:
           if fields[3] == bh_step:
              outfile_name = str(step)+"_"+str(node)+"_"+str(bh_step)+".xyz"
              outfile = open(outfile_name, 'w')
              outfile.write("%d\n" %(atom_numb))
              outfile.write("   \n")
              for i in range(atom_numb):
                  outfile.write(f.readline()) 
              outfile.close()
              print "geometry is stored in file:",outfile_name 
              return True
    f.close()

def read_exafs(filename, step, node, bh_step, spec_numb):
    f = open(filename, 'r')
    count = 0
    while True:
        line = f.readline()
        if not line:
           break
        fields = [ field for field in line.split()]
        if 'step:' in line:
           if fields[1] == bh_step:
              count+=1
              outfile_name = str(step)+"_"+str(node)+"_"+str(bh_step)+fields[3]+".dat"
              outfile = open(outfile_name, 'w')
              line = f.readline()
              fields = [ field for field in line.split()]
              while len(fields)==2:
                 outfile.write(line)
                 line = f.readline()
                 fields = [ field for field in line.split()]
        if count == spec_numb:
           f.close()
           outfile.close()
           return True


def main():
    arg = sys.argv
    #inputfile = arg[1]
    energy = arg[1]
    s = arg[2]
    if len(arg) < 4:
       steps=20
       nnode=5
       spec_numb=1
    else:
       steps=int(arg[3])
       nnode = int(arg[4])
       spec_numb = int(arg[5])
    if len(arg) < 6:
       spec_numb = 1
    else:
       spec_numb = int(arg[5])
    
    for step in range(steps):
        atoms_found = False
        for node in range (nnode):
            pot_file = "pot_log"+"_"+str(step)+"_"+str(node)
            print 'Search in ', pot_file
            bh_step = read_dots(pot_file, energy, s)
            if bh_step:
               print step, node, bh_step
               sys.exit()
            #if bh_step is not None:
            #   traj_file = "/ssd/users/leili/"+"_"+str(step)+"_"+str(node)
            #   atoms_found = read_atoms(traj_file, step, node, bh_step)
            #   exafs_file = 'exafs_'+str(step)+'_'+str(node)
            #   atoms_found = read_exafs(exafs_file, step, node, bh_step, spec_numb)
            #   break
        #if atoms_found:
        #   break

    print "not found"
if __name__ == '__main__':
    main()

        
