import numpy as np
import matplotlib.pyplot as plt
import subprocess
from data import data


##################################################################
def make_spec(fo,fsamp,amp,w):

##A function to return a spectrum of densities at the regulalrly spaced
##vector of frequencies fo due to the sampling of mode frequencies fsapm
##weighted by the vector amp. The sopectrum is created by representign
##each sample frequency as a lorentzian of width w, and summing the
##lorentzians for all samples.
##
##fo and fsamp must both be row vectors

    #print  "fo shape: ", np.shape(fo)
    #print "fsamp shape: ",  np.shape(fsamp)

    n = len(fsamp)
    if n == 1:
        fsamp = np.ones(np.shape(fsamp))
    #print "n: ", n
    z = np.shape(fo)

    s = np.zeros(z)
    i = 0 
    while True:
        r = np.subtract(fo,fsamp[i])**2
        #print "r shape: ", np.shape(r)
        s = s + amp[i]/((r + w)*np.pi*w)
        i += 1
        if i == n:
            break

    return s
##
##samp * fb = fsamp * fo;
##eb = e * ones(size(fo));
##
##spec=sum(eb*sqrt(br)./(pi*(1+*fb)),1);


##################################################################

def weight_spec(fo,s,T):

    h = 4.136e-15       # eV.s
    Kb = 8.617e-5       # eV/K
    conv = h*1e12/(Kb*T)
    BEf = h*fo/(np.exp(conv*fo)-1)
    ws = s*BEf
    return ws

############################################################################
def round_data(infile):
    
    dt = data(infile)
    outfile = infile + '-single-prescission'
#print rs   
    dt.headers["xlo xhi"] = np.round(dt.headers["xlo xhi"],5)
    dt.headers["ylo yhi"] = np.round(dt.headers["ylo yhi"],5)

    dt.headers["zlo zhi"] = np.round(dt.headers["zlo zhi"],5)
    #dt.headers["xy xz yz"] = np.round(dt.headers["xy xz yz"],6)
    x = np.round(np.array(dt.get("Atoms")),5)
    dt.replace("Atoms", 3, x[:,2])
    dt.replace("Atoms", 4, x[:,3])
    dt.replace("Atoms", 5, x[:,4])

    x = np.round(np.array(dt.get("Velocities")),5)
    dt.replace("Velocities", 2, x[:,1])
    dt.replace("Velocities", 3, x[:,2])
    dt.replace("Velocities", 4, x[:,3])

    dt.write(outfile)

############################################################################
def dm(k,Rv,Kv):
    Kv1 = np.copy(Kv)
#
# A function that returnes the dynamical matrix for wave vector k build
# using the list of stiffness matrices Kv at positions Rv
#
#
    nK, nm, nr = np.shape(Kv1)
    wv = np.exp(np.dot(Rv,k)*np.sqrt(-1 + 0j))

    Dk = np.zeros(nm)
    
    i = 0
    while True:
        Dk = Dk - np.multiply(wv[i],np.squeeze(Kv1[i]))
        i += 1
        if i == nK:
            break
    return Dk

#############################################################################
def compute_stiffness_matrix():

## Compute stiffness data

# run lammps to compute forces on each atom due to a displacment from 
# every other atom in a system of 3x3x3 unitcells of Si
#!./lmp < in.Si-hessian
    """# Run LAMMPS and compute the energy of the strined Cu sample
    p = subprocess.Popen("lmp < in.minimize-eam", stdout=subprocess.PIPE, shell=True)       #./lmp < in.minimize-eam   
    (output, err) = p.communicate()
    #print output """

# Convert the output from lammps into a global mass-wieghted stiffnes matrix 
# (Hessian matrix) 
#Run 
    p = subprocess.Popen("./forcedump2hessian.py dump.forces hessian.d", stdout=subprocess.PIPE, shell=True) 
    (output, err) = p.communicate()
#print output 

# Start by loading the data required to build the dynamical matrix for Si. 
# This consists of 27 24x24 stifness matrixes and a list of position 
# vectors. They are stored in the directory Stiffness_Data

# Lattice parameter
    a = 5.43 # Ang

    na = 8  # Number of atoms per unit cell
    nm = 24 # Number of degrees of freedom per unit cell
    nK = 27 # number of stiffness matreses that contribute to the (3x3x3=27)
         # dynamical matrix         ln ~
    Rv = []
# Define the position vectors of the corner of each of the 27 unit cells
    """for cell = 1:nK
        [i,j,k] = ind2sub([3,3,3],cell)
        Rv(cell,:) =  ([i,j,k]-2)*a
    end     """

    cell = 0
    while True:
        ar = np.array(np.unravel_index(cell,(3,3,3)))
        Rv.append(list((ar-2)*a))          #create Rv np.unravel converts index to i,j 
        cell += 1                                                      #ex. 5 in a 3x3 is (1,2)
        if cell == nK:
            break
    #print "shape Rv", np.shape(Rv)


                                    # cells from the central unit cell
    kmax = np.pi/a                    # Distance in reciprical space to the 
                                    # Brillouin zone edge   
    h = []
# Construct the list of stifness matricies
# Load the Hessian for the full 27 unit cell system
    with open ('hessian.d', "r") as myfile:
        for line in myfile:
            h.append(list(np.fromstring(line, dtype=float, sep=' ')))
    #print type(h)
    #print np.shape(h)
    h = np.array(h)
    h = 0.5 * (h + h.transpose())

# cut out the central column of sub matricie
    rngj = [(nK//2)*24,(nK//2+1)*24]
    #print "rngj:", rngj
    Kv = np.zeros((nK, nm, nm))
    i = 0
    while True:
        rngi = [(i)*24,(i+1)*24]
        #print rngj
        #print "rngi: ", rngi
        #print h[rngi[0]:rngi[1]]
        #print h[rngi][rngj]
        #print h[rngi,rngj]
        hT = h[rngi[0]:rngi[1],rngj[0]:rngj[1]] #np.transpose(
        #print "hT shape", np.shape(hT)
        #print hT.tolist()
        Kv[i] = hT.tolist()
        #print Kv[i]
        #print "kv shape:", np.shape(Kv)
        i += 1                                                      #ex. 5 in a 3x3 is (1,2)
        if i == nK:
            break      
    '''# Construct the list of stifness matricies
    i   i = 0
    f   while True:
        s = 'Stiffness_Data/K_' + str(np.ravel(i))          #num2str(i));
        Kv(i,:,:) = load(s);
        i += 1                                                      #ex. 5 in a 3x3 is (1,2)
        if i == nK:
            break'''


    # Test that these work: compute the frequencies at the Gamma point (the
    # point where wavelength is infinite, i.e. block deformation
    k = np.array([0,0,0])
    Rv = np.array(Rv)
    Kv = np.array(Kv)
    Dk = dm(k,Rv,Kv)   # input np arrays
    #print Dk
    #matrix stuff for ploting
    Dk_plt = np.copy(Dk)
    nr, nc = Dk_plt.shape 
    extent = [-0.5, nc-0.5, nr-0.5, -0.5] 

    # Plot Dk
    ax = plt.subplot2grid((1,1),(0,0))
    ax.imshow(np.real(Dk_plt), extent=extent, origin='upper', interpolation='nearest')
    #ax.set_ylabel('i')
    #ax.set_title(Title)

    wg = wsort(Dk)
    #print "wg= ", wg

    # Clear old data
    #clear('wtab','ktab','kpath');
    #plt.show()
    return k, Rv, Kv, Dk, wg, a, nm

def wsort(mat):
    #print "mat:", mat
    mat1 = np.copy(mat)
    mat1 = mat1 + 0j #Dk1[Dk1<0] = 0
    w = np.sort(np.real(np.sqrt(mat1)))
    return w

################################################################
def heaviside(x):
    return 0.5 * (np.sign(x) + 1)

################################################################

# if this script is run as standalone do stuff
if __name__ == "__main__":
    infile  = "./data.Si-equilibrated-300K"
    round_data(infile)
    #compute_stiffness_matrix()
