import sys
import os
from ase.io import read, write

def main():
    arg = sys.argv
    inputfile = arg[1]
    fmt1 = arg[2]
    outputfile = arg[3]
    fmt2 = arg[4]
    p1 = read(filename=inputfile, index=0, format = fmt1)
    p1.center()
    for atom in p1:
        if atom.position[2] < 2.0:
           print("%s   %d" % ("  LIST", atom.index+1))
    write(filename=outputfile, images=p1, format = fmt2)

if __name__ == '__main__':
    main()

