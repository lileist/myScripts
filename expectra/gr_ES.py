#!/usr/bin/env python
"""
This code is used to analyze the distribution of energy and s (chi_deviation).
How to use:
    gr_ES [inputfile_name] [de] [ds]
"""
import sys
import os

def read_dots(filename):
    f = open(filename)
    energy = []
    s =[]
    accept=[]
    for line in f:
        if line.startswith('#'):
           continue
        fields = [ field for field in line.split()]
        energy.append(float(fields[6]))
        s.append(float(fields[7]))
        accept.append(float(fields[2]))
        #dot = [float(fields[6]), float(fields[7])]
    f.close()
    return energy, s, accept

def main():
    arg = sys.argv
    inputfile = arg[1]
    de = float(arg[2])
    ds = float(arg[3])
    outfile = arg[4]
    energy,s,accept = read_dots(inputfile)
    e_min = min(energy)
    e_max = max(energy)
    s_min = min(s)
    s_max = max(s)
    e_gr=[]
    s_gr=[]
    e_outfile= 'e_gr'+'_'+outfile+'.dat'
    s_outfile= 's_gr'+'_'+outfile+'.dat'
    log_e = open(e_outfile,'w')
    log_s = open(s_outfile, 'w')

    #initialize e_gr and s_gr
    for i in xrange (int((e_max - e_min)/de)+1):
        e_gr.append(0.0)
    for i in xrange (int((s_max - s_min)/ds)+1):
        s_gr.append(0.0)

    #analyze
    numb_accept = 0.0
    for i in xrange (len(energy)):
        numb_accept = numb_accept + accept[i]
        for j in xrange (int((e_max - e_min)/de)+1):
            if energy[i] >= e_min + float(j)*de and energy[i] < e_min + float((j+1))*de:
               e_gr[j] += 1
        for k in xrange (int((s_max - s_min)/ds)+1):
            if s[i] >= s_min + float(k)*ds and s[i] < s_min+float((k+1))*ds:
               s_gr[k] += 1
    print "acceptace: ", numb_accept / float(len(accept))
    #normalize
    for i in xrange (int((e_max - e_min)/de)+1):
        e_gr[i] = float(e_gr[i])/float(len(energy))
        r = e_min + float(i)*de - de*0.5
        log_e.write('%15.6f  %15.6f\n' % (r, e_gr[i]))
    for i in xrange (int((s_max - s_min)/ds)+1):
        s_gr[i] = float(s_gr[i])/float(len(s))
        r = s_min + float(i)*ds - ds*0.5
        log_s.write('%15.6f  %15.6f\n' % (r, s_gr[i]))
    log_e.close()
    log_s.close()
if __name__ == '__main__':
    main()

        
