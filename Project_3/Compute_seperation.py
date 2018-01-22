import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def get_frame(filename):
    #num_lines = sum(1 for line in open('myfile.txt'))    #get file lenth (not needed) 
    with open(filename, 'r') as f:
        print "running get"
        j = []
        i = 0
        frames = []
        for line in f:
            i+=1
            if "ITEM: ATOMS" in line:
                j.append(i)

        for line in j:
                print "testing"
                frames.append(np.genfromtxt(filename, skip_header=(line), skip_footer=(i-(line+216))))
                print np.shape(frames)
        return frames

def phasespace_seperation():
    "runing phase"
    dt = 1e-2       #dt used in lamps in picoseconds (1e-3 = 1 femptosecond)
    #ht = 50        #heating time (simulation time in picoseconds)
    sf = 1000       #sampeling frequency
    filenameA = 'pv-FS.d'
    filenameB = 'pv-FD.d'
    sA = get_frame(filenameA)
    print "sA: ", np.shape(sA)
    sB = get_frame(filenameB)
    print "sB: ", np.shape(sB)

    i = 1
    FrLen= len(sA)
    d = np.zeros((FrLen,1))
    a = np.zeros((FrLen,1))
    while True:
        sA1=np.array(sA[i])
        sB1=np.array(sB[i])
        #print "sA-shape: ", np.shape(sA1)        
        sxA = np.hstack((sA1[:,1],sA1[:,2],sA1[:,3]))
        #sxA = sA1[:,1:4]
        #print "sxA ", sxA
        sxB = np.hstack((sB1[:,1],sB1[:,2],sB1[:,3]))
        #sxB = sB1[:,1:4]
        #print "sxB ", sxB
        v = sxA - sxB
        print "v-shape: ", np.shape(v)
        #print "v: ", v
        
        ddot=np.dot(v,v.T)
        print ddot
        d[i] = np.sqrt(ddot)
        print "d[i]: ", d[i]   

        svA = np.hstack((sA1[:,4],sA1[:,5],sA1[:,6]))
        svB = np.hstack((sB1[:,4],sB1[:,5],sB1[:,6]))
        #svA = sA1[:,4:7]
        #svB = sB1[:,4:7]

        a[i] = np.dot(svA,svB)/np.sqrt(np.dot(svA,svA)*np.dot(svB,svB))
        print "a[i]: ", a[i]
        i +=1
        if i >= FrLen:
            break

    a = np.arccos(a)*180/np.pi
    x = np.arange(0,FrLen)*sf*dt

#calc Lyapunov instability

    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    nand = np.nan_to_num(d)    
    nand[nand == 0] = 10**-10
    #nand = np.log10(d)
    nand = nand[:,0]
    print "nand_size: ", np.shape(nand)
    popt, pcov = curve_fit(func, x, nand )
    print "popt: ", popt
    print "pcov: ", pcov
#plot curves
    ax = plt.subplot2grid((1,2),(0,0))
    ax.plot(x, np.log10(nand), label="Simulation Data")
    ax.plot(x, func(x, *popt), label="Curve fit")
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Log[Phasespace seperation (Ang)]')
    ax.legend()

    ax1 = plt.subplot2grid((1,2),(0,1))
    ax1.plot(x, a)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Angle between trajectories (deg)')

    plt.show()
#############################



if __name__ == "__main__":
    phasespace_seperation()
    #filename = 'pv.d'
    #frames  = get_frame(filename)
    #print "frams[0]: ", frames[0]
    #print "frames[1]: ", frames[1]
