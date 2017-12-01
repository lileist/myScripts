#!/usr/bin/env python
import sys
from expectra.atoms_operator import match, single_atom
from ase.io import read
def main():
    arg = sys.argv
    p1=read(arg[1],format='xyz')
    p2=read(arg[2],format='xyz')
    matched=match(p1, p2, 0.01, 3.0, indistinguishable=True, )
    if matched:
       print "matched"
    print matched
if __name__ == '__main__':
    main()
