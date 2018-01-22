import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack 
from compute_tools import *

def CreateSpec(Temp):
    file = 'vacf.d-' + str(Temp) + 'K'

# Load the VACF data
    m = np.loadtxt(file)

    # scale time step
    dt = 1.0e-3
    t = m[:,0]*dt/2     #not sure whay it is needed bu the /2 make it plot right
    fs = 24

    #x = np.arange(0,FrLen)*sf*dt

    ax = plt.subplot2grid((1,2),(0,0))
    ax.plot(t,m[:,3],'b-')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('VACF (Arb. U)')


    npi = len(m)
    print "npi= ", npi
    #np = 100;
    fmax = 40

    # The manual (slow) way to do it
    f = []
    g = []
    i = 0
    """while True:
        f.append(i*fmax/npi)
        c = np.cos(f[i]*2*np.pi*t)
        g.append(np.trapz(c*m[:,3])**2)
        i += 1 
        print "i= ", i
        if i == npi :
            break"""
    # the fast way to do it 
    f   = np.arange(0,(npi-1))/(2*(t[1]-t[0])*(npi-1))
    g   = fftpack.dct(m[:,3], norm='ortho')**2
    print "g-shape: ", np.shape(g)

    rng = np.nonzero(heaviside(fmax-f))
    print "rngs: ", np.shape(rng)
    # Build a spectrum
    fo = np.arange(1, fmax, 0.01)
    w  = 0.01   
    s  = make_spec(fo,f[rng],g[rng],w)
    print np.shape(g[rng])
    s  = s/np.max(s)

    #h = plot(f(rng),g(rng)/max(g(rng)),'b-',fo,s,'r-')
    ax1 = plt.subplot2grid((1,2),(0,1))
    ax1.plot(fo,s,'r-')
    ax1.set_xlabel('Frequency (THz)')
    ax1.set_ylabel('Intensity (Arb. U)')

    plt.show()

CreateSpec(900)
