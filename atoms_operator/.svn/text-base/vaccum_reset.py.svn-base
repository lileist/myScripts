#!/usr/bin/env python
import sys
import os
from ase.io import read, write
from ase.io.vasp import write_vasp


def main():
    arg = sys.argv
    inputfile = arg[1]
    fmt1 = arg[2]
    c = float(arg[3])
    a = None
    b = None
    if len(arg)==6:
       a = float(arg[4])
       b = float(arg[5])
    outputfile = 'POSCAR'
    p1 = read(filename=inputfile, index=0, format = fmt1)
    if a is None and b is None:
       a = p1.get_cell()[0]
       b = p1.get_cell()[1]
       p1.set_cell([a ,b ,[0.,0., c ]],scale_atoms=False)
    else:
       p1.set_cell([[a,0.0,0.0] ,[0.0,b,0.0] ,[0.,0., c ]],scale_atoms=False)
    p1.center()     
    write_vasp(filename = outputfile, atoms = p1, direct=True, vasp5=True)
if __name__ == '__main__':
    main()
