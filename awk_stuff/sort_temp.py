import fileinput
import os

for line in fileinput.input("Example.cif.dat.2nd.txt", inplace=1):
    line = line.rstrip(os.linesep)
    print(" ".join(sorted(line.split())))
