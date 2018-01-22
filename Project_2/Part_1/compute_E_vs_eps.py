# A Python script to compute strain energy as a function of  strain
# amplitude et
import subprocess
from strain_tensor import strain_tensor 
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from LaTexConv import poly2latex
from view_strain import view_strain
av = []
ev = []

eta = -0.05       #start 
step = 0.001
stop = 0.05
#mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])          #hydrostatic strain
#mat = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])          #C11 - Uniaxial strain - z direction (e33)
#mat = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]])          #C12 - Biaxial strain - z,y (e22e33)
mat = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]])           #C44 - Shear strain - (e12)
#mat = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]])         # Example of a pure shear deformation in the x direction on the y plane (e11(-e22))
view_strain(mat*0.1)        #exagerated for easy viewing
while True:
    av.append(eta)

    strain_tensor(eta * mat) 
    # Run LAMMPS and compute the energy of the strined Cu sample
    p = subprocess.Popen("lmp < in.minimize-eam", stdout=subprocess.PIPE, shell=True)       #./lmp < in.minimize-eam   
    (output, err) = p.communicate()
    #print output 

    # Extract the data from the LAMMPS output file log.lammps
    p = subprocess.Popen("./get_energy.sh", stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    #print "NRG=", output

    # Read that value in to the energy vector
    ev.append(float(output))            #convert to float from string
    
    #step the while
    eta += step

    #break while if stop it reached
    if eta > stop:
        break
#convert ev, av to nparrays  - prevent errors
av = np.asarray(av)
ev = np.asarray(ev)

# Use the polyfit function to return the coefficients of the 2nd,
# 3rd and 4th order polynomials that best fit the energy data.
c2 = np.polyfit(av,ev,2)
p2 = np.poly1d(c2)          #for convienince
f2 = p2(av)                 #create aray of fit points
#print p2
R2 = stats.linregress(ev, f2)        #get R squared
#print "%.8f" % R2[0]

c3 = np.polyfit(av,ev,3)
p3 = np.poly1d(c3)
f3 = p3(av)                 #create aray of fit points
#print str(p3)
R3 = stats.linregress(ev, f3)        #get R squared
#print "%.8f" % R3[0]

c4 = np.polyfit(av,ev,4)
p4 = np.poly1d(c4)
f4 = p4(av)                 #create aray of fit points
#print str(p4)
R4 = stats.linregress(ev, f4)        #get R squared
#print "%.10f" % R4[0]


#Plot It all up
ax = plt.subplot2grid((1,1),(0,0)) 
ax.plot(av, ev, 'ro-', label = "Simulation Data") 
ax.plot(av, f2, 'b-', label = poly2latex(c2, width = 5) + "\n$R^2 =  %.8f$" % R2[0]) 
ax.plot(av, f3, 'k-', label = poly2latex(c3, width = 5) + "\n$R^2 =  %.8f$" % R3[0]) 
ax.plot(av, f4, 'g-', label = poly2latex(c4, width = 5) + "\n$R^2 =  %.10f$" % R4[0]) 
ax.set_xlabel('$\eta$')
ax.set_ylabel('Energy')
ax.legend(loc=9)    
#plt.tight_layout()
plt.show()

# The energy units from LAMMPS are in eV, and distances are given in
# Angstroms. We need to convert the energies to Si to compute
# elastic constants in GPa 
eV = 1.60217646e-19     #[eV]
conv = eV*1e30          #[eV*m
# lattice parameter of the Cu unit cell (when  eta = 0) 
ao = 3.615


# Using the second order polynomoial fit:
# Compute the value for eta the gives the minimum energy
p2_prime = np.polyder(p2)
xo = fsolve(p2_prime,-0.05)
print xo
a = (1+xo) * ao    # relaxed lattic paramter
print "relaxed lattice (2nd order) = ", a
# Depending on how we are straining the copper unit cell we can now 
# relate the second derivative of the energy to the elastic constants.
# For the case of a uniaxial tension we can relate C_11 to the energy 
k2 = p2_prime(0)  # Compute Stiffness
p2_2prime = np.polyder(p2_prime)  
k2 = p2_2prime(0) # Compute Stiffness
print "stiffness (2nd order)", k2
k2 = conv*k2/(a**3) # Convert to units of Pa
print "stiffness (2nd order) %s [GPa]" %(k2/(10**9))

# Using the third order polynomoial fit:
# Compute the value for eta the gives the minimum energy
p3_prime = np.polyder(p3)
xo = fsolve(p3_prime,-0.05)
print xo
a = (1+xo) * ao    # relaxed lattic paramter
print "relaxed lattice (3nd order) = ", a
# Depending on how we are straining the copper unit cell we can now 
# relate the second derivative of the energy to the elastic constants.
# For the case of a uniaxial tension we can relate C_11 to the energy 
k3 = p3_prime(0)  # Compute Stiffness
p3_2prime = np.polyder(p3_prime)  
k3 = p3_2prime(0) # Compute Stiffness
print "stiffness (3nd order)", k3
k3 = conv*k3/(a**3) # Convert to units of Pa
print "stiffness (3nd order) %s [GPa]" %(k3/(10**9))

# Using the fourth order polynomoial fit:
# Compute the value for eta the gives the minimum energy
p4_prime = np.polyder(p4)
xo = fsolve(p4_prime,-0.05)
print xo
a = (1+xo) * ao    # relaxed lattic paramter
print "relaxed lattice (4nd order) = ", a
# Depending on how we are straining the copper unit cell we can now 
# relate the second derivative of the energy to the elastic constants.
# For the case of a uniaxial tension we can relate C_11 to the energy 
p4_2prime = np.polyder(p4_prime)  
k4 = p4_2prime(0) # Compute Stiffness
print "stiffness (4nd order)", k4
k4 = conv*k4/(a**3) # Convert to units of Pa
print "stiffness (4nd order) %s [GPa]" %(k4/(10**9))


