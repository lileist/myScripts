#!/usr/bin/env python
"""
This code is used to extract geometry and exafs for given energy and s (chi_deviation).
How to use:
   dot_selection  [energy] [s] [pl_steps] [bh_steps] [species_number]
"""
import sys
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_pl(filename):
    f = open(filename, 'r')
    pl_cycle = None
    node_numb=0
    pl =[]
    dots=[]
    while True:
        line = f.readline()
        fields = line.split()
        if not line:
           break
        if "======" in line:
           if pl_cycle is not None and len(dots) >0:
              pl.append(copy.deepcopy(dots))
           continue

        if "pl_cycle" in line:
           pl_cycle = int(fields[1].replace(",",""))
           node_numb = int(fields[3])
           continue

        #if "replace" or "insert" or "pop" or "append" in line:
        if ":" in line:
           #new dots found. paretoLine needs to be updated
           dots=[]
           continue
        if len(fields) == 3:
           dots.append([float(fields[1]),float(fields[2])])
    f.close()
    return pl

def read_vc(filename):
    f = open(filename, 'r')
    configs=[]
    visited_configs=[]
    while True:
        line = f.readline()
        fields = line.split()
        if 'pareto_step' in line:
           pl_step = int(fields[1])
           configs.append(visited_configs)
           visited_configs=[]
           continue
        if not line:
           break
        #visited_configs[ fields[0] ]= [float(fields[1]), float(fields[2])]
        visited_configs.append( [float(fields[1]), float(fields[2])] )
    return configs

def read_accepted_dots(filename, pl_step, nnode):
    accepted_dots =[]
    dots=[]
    for i in range(pl_step):
        #dots=[]
        for j in range(nnode):
            pot_file = filename+'_'+str(i)+'_'+str(j)
            f = open(pot_file, 'r')
            while True:
                line=f.readline()
                if not line:
                   break
                if "#" in line:
                   continue
                fields=line.split()
                if fields[2]=='1' and fields[1]!='-1':
                   dots.append([float(fields[6]),float(fields[7])])
        if len(dots)>0:
            accepted_dots.append(copy.deepcopy(dots))
    return accepted_dots

def extract_xy(data_set):
    x_set=[]
    y_set =[]
    for dots in data_set:
       x_data =[]
       y_data=[]
       for dot in dots:
           x_data.append(dot[0])
           y_data.append(dot[1])
       x_set.append(np.array(x_data))
       y_set.append(np.array(y_data))
    return np.array(x_set), np.array(y_set)
    
args=sys.argv
if len(args)>1:
   vc_dots = read_accepted_dots(os.getcwd()+'/'+'pot/pot_log', int(args[1]), int(args[2]))
else:
   vc_dots = read_vc("visited_configs.dat")
print len(vc_dots)
pl_dots=read_pl("paretoLine.dat")

x1_data_set, y1_data_set=extract_xy(pl_dots)

x_data_set, y_data_set=extract_xy(vc_dots)

fig, ax = plt.subplots()
#initialize of line object
line, = ax.plot(np.random.rand(100), 'bo')
line2, = ax.plot(np.random.rand(100), color='red',linestyle='solid', linewidth=2)
ax.set_xlim(-188, -185)
#ax.set_aspect(aspect=0.25, adjustable='box')
#ax.set_aspect(aspect='equal')
ax.set_ylim(11, 23)
def update(i):
    #update line object with data
    line.set_ydata(y_data_set[i])
    line.set_xdata(x_data_set[i])
    line2.set_ydata(y1_data_set[i])
    line2.set_xdata(x1_data_set[i])
    return [line, line2]

plt.plot([-188,-185],[12.17913,12.17913],"k--",lw=2)
plt.plot([-187.891137,-187.891137],[11,23],"k--",lw=2)
ani = animation.FuncAnimation(fig, update, np.arange(0, len(x_data_set)), interval=100)
#animation.save(filename='movie')
plt.show()


