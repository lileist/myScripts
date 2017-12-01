#!/usr/bin/env python
"""
This code is used to extract geometry and exafs for given energy and s (chi_deviation).
How to use:
   dot_selection  [energy] [s] [pl_steps] [bh_steps] [species_number]
"""
import sys
import os

def read_dots(filename):
    f = open(filename)
    dot = []
    for line in f:
        if line.startswith('#'):
           continue
        fields = [ field for field in line.split()]
        #print fields[6], " ", fields[7]
        dot.append([fields[6],fields[7], fields[5]])
    f.close()
    return dot

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
    steps=int(arg[1])
    nnode = int(arg[2])
    dots = []
    for step in range(steps):
        for node in range (nnode):
            pot_file = "pot_log"+"_"+str(step)+"_"+str(node)
            if os.path.isfile(pot_file):
               dots.extend(read_dots(pot_file))
            else:
               continue
    log_dots = open('t_pot.dat', 'w')
    log_alpha = open('alpha.dat','w')

    alpha_list=[]
    prob = []
    interval = 1.0 / float(nnode)
    for i in range (nnode):
        alpha_list.append(float(i) * interval)
        prob.append(0)
    alpha_list.append(1.0)

    alpha = float(dots[0][2])
    alpha_sum = alpha
    alpha_avg = []
    index = 1
    prob[0] = 1
    for i in range(len(dots)):
        if float(dots[i][2]) != alpha:
           alpha = float(dots[i][2])
           alpha_sum = alpha_sum + alpha
           index += 1
           if index > 50:
              for j in range(len(alpha_list)-1):
                 if alpha >= alpha_list[j] and alpha <= alpha_list[j+1]:
                    prob[j] += 1
        if index % 5 == 0 and alpha_sum != 0.0:
           alpha_avg.append((alpha_sum/5.0))
           alpha_sum = 0.0
           
        log_dots.write("%d  %15.6f %15.6f %8.4f \n" % (i, float(dots[i][0]), float(dots[i][1]), float(dots[i][2])))
    for i in range(len(alpha_avg)):
        log_alpha.write("%d %8.4f \n" % (i, alpha_avg[i]))
    for i in range(len(prob)):
        print alpha_list[i+1], float(prob[i])/ float(index-50)
    log_dots.close()

if __name__ == '__main__':
    main()

        
