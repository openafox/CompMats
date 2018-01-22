import numpy as np
import matplotlib.pyplot as plt
import subprocess

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

#compute_stiffness_matrix()
