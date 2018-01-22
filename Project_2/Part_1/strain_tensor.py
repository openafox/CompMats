#!/usr/bin/python


import sys,os
import math as m
import operator as o
import matplotlib.pyplot as plt
import numpy as np
from view_strain import *
from data import data
#print rs
def strain_tensor(s):    

# Load a data file
    dt = data("./data.Cu-unit");

    (xlo, xhi) = dt.headers["xlo xhi"]
    (ylo, yhi) = dt.headers["ylo yhi"]
    (zlo, zhi) = dt.headers["zlo zhi"]
    (xy, xz, yz) = dt.headers["xy xz yz"]
# strain the box boundaries
    #print xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    org = np.array([xlo,       ylo,     zlo      ])
    #print "org: %s" % org
    av  = np.array([xhi-xlo,   0.0,     0.0      ])
    bv  = np.array([xy,        yhi-ylo, 0.0      ])
    cv  = np.array([xz,        yz,        zhi-zlo])

    org = org + np.dot(org,s)
    #print "org: %s" %org

    av  = av  + np.dot(av,s)
    bv  = bv  + np.dot(bv,s)
    cv  = cv  + np.dot(cv,s)

# First rotate so that the "av" vector is parallel to (1,0,0) 
    ax1   = np.cross(np.array([1.0,0.0,0.0]),av)
    d     =  m.sqrt(np.dot(av,av))
    A     =  m.sqrt(np.dot(ax1,ax1))
    ang1  = -m.asin(A/d)

    #print 'ax1=',ax1
    #print 'ang1=',ang1

    if (np.dot(ax1,ax1) != 0.0):
        org  = qrot(org, ax1, ang1)
        av   = qrot(av,  ax1, ang1)
        bv   = qrot(bv,  ax1, ang1)
        cv   = qrot(cv,  ax1, ang1)
    
    #print "bv size", np.shape(bv)
    #print 'av=',av
    #print 'bv=',bv
    #print 'cv=',cv

# Now rotate about (1,0,0) so that bv vector is in x-y plane 
    #v2   = np.multiply(np.array([0.0, 1.0, 1.0]), bv)

    v2   = np.multiply(np.array([[0.0, 1.0, 1.0]]), bv)
    #print "v3", v3
    #print np.shape(v3)
    #print "v2", v2
   #print np.shape(v2)

    v2  = v2/np.sqrt(np.einsum('jk, jk', v2, v2)) #
    #print "v3", v3
    #print np.shape(v3)

    #v2 = v2/m.sqrt(np.dot(v2,v2))
    #print "v2", v2
    #print np.shape(v2)
    ang2 = -m.acos(np.dot(v2,np.array([0.0, 1.0, 0.0])))
    ax2  = np.array([1.0, 0.0, 0.0])
    #print 'v2 =',v2
    #print 'ax2=',ax2
    #print 'ang2=',ang2

    org  = qrot(org, ax2, ang2)
    #print 'org=',org
    av   = qrot(av,  ax2, ang2)
    #print 'av=',av
    bv   = qrot(bv,  ax2, ang2)
    #print 'bv=',bv  
    cv   = qrot(cv,  ax2, ang2)
    #print 'cv=',cv  



    xlo = org[0,0] 
    ylo = org[0,1] 
    zlo = org[0,2]
    xhi = av[0,0] + org[0,0]
    yhi = bv[0,1] + org[0,1] 
    zhi = cv[0,2] + org[0,2]
    xy  = bv[0,0]
    xz  = cv[0,0]
    yz = cv[0,1]

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

    dt.write("./data.Cu-unit-stretched")

#a = 0.3
#mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#mat = np.array([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
#s = a * mat 
#view_strain(s)
#print s
#strain_tensor(s)

