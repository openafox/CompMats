import numpy as np
import matplotlib.pyplot as plt
import subprocess
from compute_tools import *

# A py-script to compute the dispersion relation along the 
# [1 0 0] direction out to the edge of the first Brillouin zone 
#

## Compute and then load the stiffness data for use in constructing the 
# dynamical matrix

k, Rv, Kv, Dk, wg, a, nm = compute_stiffness_matrix()

def ToureBZ(stp, its, k, kstep, ktab, kdir, kpath, wtab, Rv, Kv):
    i = 0
    while True:
        k = ktab[stp] + kstep*kdir                      # Define wavevector k
        #print "K:", k
        ktab.append(k)      #(stp,:) = k(:)
        kpath.append(kpath[stp] + kstep)  #(stp)  = kpath(stp-1) + kstep
        Dk = dm(k,Rv,Kv)                                    # Construct dynamical matrix of k
        #print "Dk: ", Dk

#don't know why absolute is needed
        w = np.sort(np.real(np.sqrt(np.absolute(np.linalg.eigvals(Dk)))))      # Compute angular frequencies
        #print "w: ", w
        #print "/pi: ", w/(2*np.pi)
        #print "wshape: ", np.shape(w)
        #w = w  * float(its)/(its-i)
        wtab.append(w/(2*np.pi))      # (stp,:) = w/(2*pi)                              # Save as real frequencies
        stp += 1
        i += 1
        if i == its - 1:
            break
    return k, ktab, kpath, wtab, stp

## Plot  the dispersion
#
# Start at the gamma point (zero wavevector) and compute the frequencies of
# vibrations with wave vectros between [000] and [100]*pi/a
#
k = np.array([0,0,0])
Dk = dm(k,Rv,Kv)
w = wsort(np.linalg.eigvals(Dk))   # Compute angular frequencies
ktab = [k]
kpath = [0]
wtab = [w/(2*np.pi)]          # Save as real frequencies
stp = 0

# Take 200 steps in k space from (0,0,0) to (1,0,0) pi/a 
its = 200
kstep = np.pi/(a*its)
kdir = np.array([1,0,0])

k, ktab, kpath, wtab, stp = ToureBZ(stp, its, k, kstep, ktab, kdir, kpath, wtab, Rv, Kv)

##Toure through other parts of the BZ###
## Take 200 steps in k space from (1,0,0)pi/a to (1,1,0) pi/a 
its = 200
kstep = np.pi/(a*its)
kdir = np.array([0,1,0])
k, ktab, kpath, wtab, stp = ToureBZ(stp, its, k, kstep, ktab, kdir, kpath, wtab, Rv, Kv)

## Take 250 steps in k space from (1,1,0)pi/a to (0,0,0) pi/a 
its = 250
kstep = np.sqrt(2)*np.pi/(a*its)
kdir = np.array([-1,-1,0])
k, ktab, kpath, wtab, stp = ToureBZ(stp, its, k, kstep, ktab, kdir, kpath, wtab, Rv, Kv)

## Take 350 steps in k space from (0,0,0)pi/a to (1,1,1) pi/a 
its = 350
kstep = np.sqrt(3)*np.pi/(a*its)
kdir = np.array([1,1,1])
k, ktab, kpath, wtab, stp = ToureBZ(stp, its, k, kstep, ktab, kdir, kpath, wtab, Rv, Kv)

## Take 200 steps in k space from (1,1,1)pi/a to (1,1,0) pi/a 
its = 200
kstep = np.pi/(a*its)
kdir = np.array([0,0,-1])
k, ktab, kpath, wtab, stp = ToureBZ(stp, its, k, kstep, ktab, kdir, kpath, wtab, Rv, Kv)

# Plot dirpersion
ax = plt.subplot2grid((1,1),(0,0))
#ax1 = plt.subplot2grid((1,2),(0,1))

i = 0
kpath = np.array(kpath)
wtab = np.array(wtab)
ktab = np.array(ktab)

#print "kpath shape: ", np.shape(kpath)
#print "wtab shape: ", np.shape(wtab[:,i])
#print "kpath: ", kpath

while True:
    #print "wtab%s: %s" %(i,wtab[:,i])
    #ax1.plot(ktab[:,1],wtab[:,i],'r-');
    ax.plot(kpath,wtab[:,i],'r-')
    i += 1
    if i == nm:
        break
#axis([-kmax,kmax,0,max(wtab(:))*1.1])
ax.set_xlabel('Path distance in reciprocal space (1/Ang)')
ax.set_ylabel('Frequency (THz)')
#plt.show()


