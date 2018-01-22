import numpy as np
import matplotlib.pyplot as plt
import subprocess
from compute_tools import *

def make_spec(fo,fsamp,amp,w):

##A function to return a spectrum of densities at the regulalrly spaced
##vector of frequencies fo due to the sampling of mode frequencies fsapm
##weighted by the vector amp. The sopectrum is created by representign
##each sample frequency as a lorentzian of width w, and summing the
##lorentzians for all samples.
##
##fo and fsamp must both be row vectors

    print  "fo shape: ", np.shape(fo)
    print "fsamp shape: ",  np.shape(fsamp)

    n = len(fsamp)
    print "n: ", n
    z = len(fo)


    s = np.zeros(z)
    i = 0 
    while True:
        r = np.subtract(fo,fsamp[i])**2
        #print "r shape: ", np.shape(r)
        s = s + 1/((r + w)*np.pi*w)
        i += 1
        if i == n:
            break

    return s
##
##samp * fb = fsamp * fo;
##eb = e * ones(size(fo));
##
##spec=sum(eb*sqrt(br)./(pi*(1+*fb)),1);

def weight_spec(fo,s,T):

    h = 4.136e-15       # eV.s
    Kb = 8.617e-5       # eV/K
    conv = h*1e12/(Kb*T)
    BEf = h*fo/(np.exp(conv*fo)-1)
    ws = s*BEf
    return ws



# A py-script to compute the density of states in the in the first Brillouin zone 
#

##Start by loading the data required to build the dynamical matrix for Si. 
##This consists of 27 24x24 stifness matrixes and a list of position 
##vectors. They are stored in the directory Stiffness_Data

#Rv = load('Stiffness_Data/R_list');  ##Position vectors to the 27 unit 
##Construct the list of stifness matricies
#for i=1:nK,
#    s = strcat('Stiffness_Data/K_',num2str(i));
#    Kv(i,:,:) = load(s);
#end



k, Rv, Kv, Dk, wg, a, nm = compute_stiffness_matrix()

                                     ##cells from the central unit cell
#a  = Rv(end,1)                       ##Lattice parameter
kmax = np.pi/a                         ##Distance in reciprical space to the 
                                     ##Brillouin zone edge   


na = 8;  ##Number of atoms per unit cell
nm = 24; ##Number of degrees of freedom per unit cell
nK = 27; ##number of stiffness matreses that contribute to the 
         ##dynamical matrix         


    ##Now we will just stor a list of frequencies as we work our way
    ##through a grid of points in 1/8th segmet of the first BZ

nts = 20
wlist = []
ktab = []
stp = kmax/nts
sv = np.arange(0,kmax + stp,stp)              ##Scale of k 
i = 0
while True:
    j = 0
    while True:
        k = 0
        while True:
            kvec = [sv[i],sv[j],sv[k]]     ##Define wavevector kvec
            ktab.append(kvec[:])
            Dk = dm(kvec,Rv,Kv)             ##Construct dynamical matrix of kvec
            w  = np.real(np.sqrt(np.absolute(np.linalg.eigvals(Dk))))      ##Compute angular frequencies
            wlist.append(np.sort(w/(2*np.pi))) ##Save as real frequencies
            k += 1
            if k == nts+1: 
                break                          
        j += 1
        if j == nts+1: 
            break
    i += 1
    if i == nts+1: 
        break

##Now compute the spectrum as a sum of lorentzians

print "wlist shape: ", np.shape(wlist)
wlist = np.ravel(wlist)
print "wlist shape: ", np.shape(wlist)

fo = np.arange(1,35,0.01)
#fo = fo.reshape(1,len(fo))
w  = 0.01
s = make_spec(fo,wlist,1,w);
s = s*8


# set up plot
ax = plt.subplot2grid((1,1),(0,0))

kb = 1.381e-23      # in J/K
h  = 6.626e-34      # in J.s

#temp Loop
T = 300
Tstep = 300
Tstop = 1000
while True:

    # For DOS wighted by BE distribution;
    #Temp = 300;
    #BEdist = 1/(np.exp(fo*1e12*h/(kb*T))-1)
    #s = s*BEdist
    s1 = weight_spec(fo, s, T)
    ax.plot(fo,s1,'-',label = "T = %sK" %T)
    T += Tstep
    if T >= Tstop:
        break

#ax.plot(fo,s,'-')
# Plot dirpersion
ax.set_ylabel('Density of states (Arb U)')
ax.set_xlabel('Frequency (THz)')
plt.legend(loc=2)
#plt.show()

######### Cv Calc ######## got this mostly from Ryan (I am not sure the units are correct!! but need to move on)
T = 5
Tstep = 50
Tstop = 1220
T_it = []
Cv =[]
while True:

    s_int = np.trapz(s)*w   
    rho = (s*24)/(a**3*s_int)
    
    print T
    np.seterr(all = 'warn') #skip over overflow error in sinh   - Avoided by starting at 5 left for info
    Cv_T = rho/(4*kb*T**2)*((h*fo*1e12)*1/np.sinh((h*fo*1e12)/(2*kb*T)))**2
    Cv.append(np.trapz(Cv_T)*w)            
    T_it.append(T)    
    T += Tstep
    if T > Tstop:
        break
ax2 = plt.subplot2grid((1,1),(0,0))
ax2.plot(T_it,Cv)
ax2.set_ylabel('Specific heat')
ax2.set_xlabel('Temperature (K)')
plt.show()



