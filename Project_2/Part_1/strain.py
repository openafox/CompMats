#!/usr/bin/python


import sys,os
import math as m
import operator as o
import matplotlib.pyplot as plt
import numpy as np

import sys
from data import data

'''if len(sys.argv) < 1:
  raise StandardError, "Syntax: strain.py <[exx eyy ezz 2*eyz 2*exz 2*exy]>"

rs = sys.argv[1]
#print rs
rs = np.fromstring(rs,dtype=float, sep=' ')
#print rs

s = np.array([[rs[0],     rs[5]*0.5,   rs[4]*0.5],\
              [rs[5]*0.5, rs[1],       rs[3]*0.5],\
              [rs[4]*0.5, rs[3]*0.5,   rs[2]]])'''
a = 0.5
mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
s = a * mat 
view_strain(s)
#print s
          
# Load a data file

dt = data("./data.Cu-unit");

(xlo, xhi) = dt.headers["xlo xhi"]
(ylo, yhi) = dt.headers["ylo yhi"]
(zlo, zhi) = dt.headers["zlo zhi"]
(xy, xz, yz) = dt.headers["xy xz yz"]

# strain the box boundaries
#print xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
org = np.array([xlo,       ylo,     zlo      ]);
av  = np.array([xhi-xlo,   0.0,     0.0      ]);
bv  = np.array([xy,        yhi-ylo, 0.0      ]);
cv  = np.array([xz,        yz,        zhi-zlo]);

org = org + np.dot(org,s)
av  = av  + np.dot(av,s)
bv  = bv  + np.dot(bv,s)
cv  = cv  + np.dot(cv,s)

#print 'av=',av
#print 'bv=',bv
#print 'cv=',cv

# Rotate coordinate system: Need av to point along [1,0,0] and bv tp be in
# the x-y plane

# First define rotation functions
# conjugate of a quaternion
def qconj (q):
    qc = np.multiply(q, [1.0, -1.0, -1.0, -1.0])
    return qc

# Multiplication of a quaternion
def qmult (q1,q2):
    q3 = np.concatenate([[q1[0]*q2[0] - np.dot(q1[1:],q2[1:])],\
                         q1[0]*q2[1:] + q2[0]*q1[1:] + np.cross(q1[1:],q2[1:])])
    #print q3
    return q3


# Quaternion multiplication
def qrot (v,ax,ang):
    
    # quaternion rotation of vectors v1 by angle a about axis
    # ax through origin

    # Normalize the axis
    ax = ax/m.sqrt(np.dot(ax,ax))
    # create quaternions
    q1 = np.concatenate([[m.cos(ang/2)], m.sin(ang/2)*ax])
    q2 = np.concatenate([[0], v])

    q3 = qmult(q1,qmult(q2,qconj(q1)))
    # extra the rotated vector
    v2 = q3[1:]
    return v2

# First rotate so that the "av" vector is parallel to (1,0,0) 
ax1   = np.cross(np.array([1.0,0.0,0.0]),av);
d     =  m.sqrt(np.dot(av,av));
A     =  m.sqrt(np.dot(ax1,ax1));
ang1  = -m.asin(A/d);

#print 'ax1=',ax1
#print 'ang1=',ang1

if (np.dot(ax1,ax1) != 0.0):
    org  = qrot(org, ax1, ang1);
    av   = qrot(av,  ax1, ang1);
    bv   = qrot(bv,  ax1, ang1);
    cv   = qrot(cv,  ax1, ang1);

#print 'av=',av
#print 'bv=',bv
#print 'cv=',cv

# Now rotate about (1,0,0) so that bv vector is in x-y plane 
v2   = np.multiply(np.array([0.0, 1.0, 1.0]), bv)
#print v2
v2   = v2/m.sqrt(np.dot(v2,v2))
ang2 = -m.acos(np.dot(v2,np.array([0.0, 1.0, 0.0])));
ax2  = np.array([1.0, 0.0, 0.0]);

#print 'v2 =',v2
#print 'ax2=',ax2
#print 'ang2=',ang2

org  = qrot(org, ax2, ang2)
av   = qrot(av,  ax2, ang2)
bv   = qrot(bv,  ax2, ang2)
cv   = qrot(cv,  ax2, ang2)

#print 'av=',av
#print 'bv=',bv
#print 'cv=',cv

xlo = org[0]; 
ylo = org[1]; 
zlo = org[2];
xhi = av[0] + org[0];
yhi = bv[1] + org[1]; 
zhi = cv[2] + org[2];
xy  = bv[0];
xz  = cv[1];
yz  = cv[2];

dt.headers["xlo xhi"] = (xlo, xhi)
dt.headers["ylo yhi"] = (ylo, yhi)
dt.headers["zlo zhi"] = (zlo, zhi)
dt.headers["xy xz yz"] = (xy, xz, yz)

# strain the atom positions
x = np.array(dt.get("Atoms"))
x[:,2:5] = x[:,2:5] + np.dot(x[:,2:5],s)

# Rotate the strained atoms into the rotated box
# First ritate about axis 1 
if np.dot(ax1,ax1)!=0:
    for ix in range(len(x)):
        x[ix,2:5] = qrot(x[ix,2:5], ax1, ang1)

    for ix in range(len(x)):
        x[ix,2:5] = qrot(x[ix,2:5], ax2, ang2)

for ix in range(2,5):
    dt.replace("Atoms", ix+1, x[:,ix])

dt.write("./data.Cu-unit-stretched");


