#!/usr/bin/python


import sys,os
import math as m
import operator as o
import matplotlib.pyplot as plt
import numpy as np
from view_strain import view_strain
import sys
from data import data           #importer of data files from lamps
          
def strain_hydrostatic(a):
#
# A function to strain the atom coordinates in a datafile along the
# uniaxial deformation path by an amount a
#

    # Load the data file

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

# Define the deformation strain (in this case a uniaxial eloongation) 
# and its amplitude
    mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    eps = a * mat 
    view_strain(eps)

# strain the box boundaries
    hn = h;
    org = [h.xlo, h.ylo, h.zlo];
    av  = [h.xhi-h.xlo, 0.0,         0.0        ]
    bv  = [h.xy,        h.yhi-h.ylo, 0.0        ]
    cv  = [h.xz,        h.yz,        h.zhi-h.zlo]

    org = org + org * eps
    av  = av  + av  * eps
    bv  = bv  + bv  * eps
    cv  = cv  + cv  * eps

    hn.xlo = org[1]
    hn.ylo = org[2]
    hn.zlo = org[3]
    hn.xhi = av[1] + org[1]
    hn.yhi = bv[2] + org[2]
    hn.zhi = cv[3] + org[3]
    hn.xy  = bv[1]
    hn.xz  = cv[1]
    hn.yz  = cv[2]


# strain the atom positions
    bn = b;
    bn.Atoms[:,2:5] = b.Atoms[:,2:5] + b.Atoms[:,2:5]*eps

    #write_data_file('data.Cu-unit-stretched',hn,bn)
    scipy.io.savemat('data.Cu-unit-stretched1',hn,bn)
#with open("Pass.csv", "w") as out_file:
        #out_file.write("%s\t%s\t%s\n" % (user, hashed, password))



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






strain_hydrostatic(0.05)

