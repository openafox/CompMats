#!/usr/bin/python

# Scane a LAMMPS output file and compute the enthalpy of each species

import sys,os
import math as m
import operator as o
import matplotlib.pyplot as plt
import numpy as np

from dump import dump

if len(sys.argv) < 2:
  raise StandardError, "Syntax requires two arguments: forcedump2hessian.py <dumpfile> <output>"

inputfile  = sys.argv[1] # Name of LAMMPS file containing position data   
outputfile = sys.argv[2]

frame   = 0

# First read one frame and then find the number of species

d = dump(inputfile, 0)     # Create dump class but don't read the file all at once

d.map(1,"id", 2,"fx", 3,"fy", 4,"fz")

hessian = [];
mass = 28.0855
delta = 0.01
eV2J =  1.60217646e-19
amu2Kg = 1e-3/6.02214199e23
conv = 2*(eV2J/amu2Kg) * (1e-4)# conversion factor from eV/(AMU)/A^2 to J//Kg/m^2 in THz^2
#conv *= 10
print 'conv = ',conv

time = d.next()
id,fx,fy,fz = d.vecs(time,"id","fx","fy","fz")
d.tselect.none()
d.delete()

nat = len(id)
nmd = nat * 3

fo = np.array([fx,fy,fz]).T

#print fo

#range = 8
for i in range(nmd):
#for i in range(1):

    time = d.next()
    id,fx,fy,fz = d.vecs(time,"id","fx","fy","fz")
    fp = np.array([fx,fy,fz]).T
    d.tselect.none()
    d.delete()
    
    time = d.next()
    id,fx,fy,fz = d.vecs(time,"id","fx","fy","fz")
    fm = np.array([fx,fy,fz]).T
    d.tselect.none()
    d.delete()

    K = conv*(fp-fm)/(2.0*delta*mass)
    hessian.append(np.concatenate(K))
#    print 'hesian = ',hessian[-1]
    
# Write the output to file
f = open(outputfile, 'w')
for i in range(len(hessian)):
    s = ' '.join(str(cell) for cell in hessian[i])
#    print s
    f.write(s+'\n')

f.close();
