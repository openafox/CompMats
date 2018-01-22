
# -*- coding: utf-8 -*-
"""
Created on Sat April 13, 2015
Author: Austin Fox

Notes: 
        1.  % = remainder,  // = floor,  ** = power
        2.  Numpy Arrays are Base 0
        3.  Numpy Arrays are shallow copied by defalt 
            ie x = y and all opperations performed on x will happen to y
            to creat a true copy you must use x = np.copy(y)
        4.  

"""

import numpy as np
import matplotlib.pyplot as plt

#Function to Wrap an Index in either direction
def wrap_i(i,n, direction=1):            #direction =1: (i + 1), -1: (i - 1)
    if direction == -1:                 # this is actually unneeded python will wrap for you in the - direction (ie mats[-1] = mats[n-1])
        i = (i + n - 1)%(n)
    else:
        i = (i+1)%(n)
    return i

#Calculate Hamiltonian (Energy)
def ham(mat,J,h):
    A=np.copy(mat)
    Arollx = np.roll(A, 1, axis=0)      #roll matrix 1 in i direction equivilent to matlab circshift
    Arolly = np.roll(A, 1, axis=1)      #roll matrix 1 in j direction equivilent to matlab circshift
    e = Arollx*A +Arolly*A
    H =-J*np.sum(np.sum(e, axis=1), axis=0) - h*np.sum(np.sum(A, axis=1), axis=0)
    return H
    
#Function to flip a spin
def flip_spin(m,i,j):
    n = np.copy(m)      #see Note 3
    n[i,j] = m[i,j]*(-1)    
    return n

#Calculate the delta Energy
def delta_E(J,h,m,i__,j__):         #this function uses index wrapping. (I think this is fastest way as you are only shifting a few indicies)
                                    #could also use 4rolls (-1 and +1 in i and -1and +1 in j)
    n = len(m)
    #initialize wrapped positions (funky formats to make function look nice below)
    ip1=wrap_i(i__,n)       #i+1
    jp1=wrap_i(j__,n)       #j+1
    in1=wrap_i(i__,n,-1)    #i-1    
    jn1=wrap_i(j__,n,-1)    #j-1    
    E=[1,1]                 #setup matrix to hold Ei,Ef
    i=0
    m2 = flip_spin(m,i__,j__)
    mega = [m,m2]
    for mat in mega:
    #create sums of all affected i.e anything that would include m[i,j]
        sum1 =    (mat[i__,j__] * mat[ip1,j__] + mat[i__,j__] * mat[i__,jp1]) \
                + (mat[in1,j__] * mat[i__,j__] + mat[in1,j__] * mat[in1,jp1]) \
                + (mat[i__,jn1] * mat[ip1,jn1] + mat[i__,jn1] * mat[i__,j__])
        sum2 = mat[i__,j__]
        E[i] = -J*sum1-h*sum2
        i+=1
    delH=E[1]-E[0]      #Ef-Ei
    return delH
    
#test delta E vs Ei-Ef from ham (Loops through all points in matrix)
def test_delta_E(J,h,mat):
    itt = 0 #initialize iterations
    k=0     #Init k
    l=0
    n = len(mat)
    while True:
        m = np.copy(mat)        #see note 3
        while True:
            print "########k: %d l: %d##########)" %(k,l) 
            Ei = ham(m, J, h)        #Calc Ei
            #print "TEi: %d" %Ei
            flpm = flip_spin(m,k,l)  #Flip the spin of (k,l)
            #print "f: %s" %flpm
            Ef = ham(flpm, J, h)        #Calc Ef
            #print "TEf: %d" %Ef    
            DE= Ef-Ei                   #Calc Delta Energy
            print "Delta E = %s" %DE
            print "Delta E func = %s" %delta_E(J, h, m, k,l)
            k +=1
            if k == n:
                k=0
                break
        l+=1
        if l == n:
            break

#Define the Monte Carlo
def MCflip_Mag(J,h,m,B,start,stop):
    it = 0
    n = len(m)
    its = []        #initalize array to contain iteration numbers
    M = []          #initalize array to contain Magnitization values
    while True:     #Loop 
        i = np.random.randint(n)        #random number for 0 to n not inclusive
        j = np.random.randint(n)
        delE = delta_E(J,h,m,i,j)       #Calc delta E
        p = np.exp(-delE*B)             #Calc probability
        r = np.random.random()          #Random number from 0 to 1
        if r < p:                       #compair r and p and keep spin if r < p
            m[i,j] = m[i,j]*-1          #could use flipspin func but actuall adds more text and will take longer
        if it > start:                  #If Itteration number is greater then the desired minimum collect iteration and matrix data
            its.append(it)              #i.e after it converges
            M.append(np.sum(np.sum(m,0),0))
        it +=1
        if it == stop:      #stop loop if reached set max
            break
    return its, M, m        #return iterations array, Magnitization aray, and final matrix


#P1_9
#Test the Convergance - Run the desired MC and plot it
def test_con(start, stop, J, h, B, mat):
    m = np.copy(mat)
    its, M, m = MCflip_Mag(J,h,m,B,start,stop)   
    
    #Plotting   (Matricies can be plotted with plt.matshow(mat) but I like it to be cleaner so extra stuff
    nr, nc = mat.shape                          #Get matrix shape
    extent = [-0.5, nc-0.5, nr-0.5, -0.5]       #use shape to calc the "extent"

    ax = plt.subplot2grid((2,2),(0,0))          #Plot it in the upper left coner
    ax.imshow(mat, extent=extent, origin='upper') 
    ax.set_xlabel('j')
    ax.set_ylabel('i')

    ax1 = plt.subplot2grid((2,2),(0,1))          #upper right
    ax1.imshow(m, extent=extent, origin='upper') 
    ax1.set_xlabel('j')
    ax1.set_ylabel('i')    
    
    ax2 = plt.subplot2grid((2,2),(1,0),colspan = 2)     #Plot across the bottom
    ax2.plot(its, M) 
    ax2.set_xlabel('iterations')
    ax2.set_ylabel('Magnitization')
    
    plt.tight_layout()
    plt.show()

#Calc analytic expression
def Mex(T):     #T is %Tc
    T1 = np.copy(T)
    T1[T1>1] = 1
    #print T
    ret = 1-np.power((np.sinh(np.log(1+np.sqrt(2))*1/T1)),-4)
    #print ret
    ret = ret.astype(complex)
    #print ret
    ret = np.power(ret,1/8.0)
    #print ret
    ret = ret.real
    return ret

#print Mex(np.array([0.5,0.7, 0.8, 0.9, 0.99, 1, 2]))

########################################
# Part 10 - 12 
def MvsT(start, stop, J, h, Beta, Tc, Ti, Tf, Tstep, n, po):
    mi = 0
#run MC 
    while True: # matrix loop
        its = []        #iterations list
        M = []          #Mag list
        mats = []       #Matricies of matricies
        T = []          #Temp list
        Ti_it = Ti        #make copy of Ti  of iteration purposes
        if mi == 0:
            mat1 = np.ones((n,n))  #create nxnarray of 1's
            mat = np.copy(mat1)
            Title = "All Up"
        elif mi == 1:
            mat1 = np.sign(np.random.random((n,n))-po)
            mat = np.copy(mat1) 
            Title = "Random"        
        
        while True:     #MC Loop
            B = Beta / (Ti_it)       #Beta/
            itsa, Ma, mat = MCflip_Mag(J,h,mat,B,start,stop) 
            M.append(Ma)
            its.append(itsa)
            mats.append(mat)
            T.append(Ti_it)

            Ti_it += Tstep
            if Ti_it >= Tf:
                break

        #matrix stuff for ploting
        nr, nc = mat.shape 
        extent = [-0.5, nc-0.5, nr-0.5, -0.5] 

        # Plot base matrix initial and final
        ax = plt.subplot2grid((3,2),(mi,0))
        ax.imshow(mat1, extent=extent, origin='upper') 
        ax.set_xlabel('j')
        ax.set_ylabel('i')
        ax.set_title(Title)

        ax = plt.subplot2grid((3,2),(mi,1))
        ax.imshow(mat, extent=extent, origin='upper') 
        ax.set_xlabel('j')
        ax.set_ylabel('i')
        ax.set_title(Title + ' final')

        if mi == 0:
            dev1 = np.array(np.std(M,axis=1))       #creats array of stddevs
            #print dev1.shape
            Mean1 = np.array(np.mean(M,axis =1))   #creates array of M means
            #print Mean1.shape
            M1 = M
            mats1 = mats
        else:
            dev2 = np.std(M,axis=1)       #creats array of stddevs
            Mean2 = np.mean(M,axis =1)    #creates array of M means


            M2 = M
            mats2 = mats
        mi+=1
        if mi == 2:
            break

    T = np.array(T)         #make T same typ as mean and dev so plotting works
    mex = Mex(T)*Mean1.max()

    #plot MvsT
    ax1 = plt.subplot2grid((3,2),(2,0), colspan = 2)
    ax1.errorbar(T, Mean1,yerr=dev1, label="All Up") 
    ax1.errorbar(T, Mean2, yerr=dev2, label='Random') 
    ax1.plot(T, mex, label = "Exact")
    ax1.set_xlabel('Temp % of Tc')
    ax1.set_ylabel('<M>')
    ax1.legend(loc=1)
    plt.tight_layout()    
    plt.show()


####### Mag susceptibility ######
    var1 = np.var(M1,axis=1) 
    var2 = np.var(M2,axis=1)
    X1=(var1/T)#.tolist() #convert to list for plotting
    X2=(var2/T)        #.tolist() #convert to list for plotting
    a1 = plt.subplot2grid((1,1),(0,0))
    a1.plot(T, X1, label='All up' ) 
    a1.plot(T, X2,label='Random' )
    a1.set_xlabel('Temp % of Tc')
    a1.set_ylabel('X')
    a1.legend(loc=1)  
    plt.show()
    
#Auto Corr
    cf1 = autoc(mats[0])
    cf2 = autoc(mats[int(len(T)*0.8)])
    cf3 = autoc(mats2[0])
    cf4 = autoc(mats2[int(len(T)*0.8)])

    b1 = plt.subplot2grid((2,2),(0,0))
    b1.matshow(cf1)
    #b1.imshow(cf1, extent=extent, origin='upper') 
    b1.set_xlabel('j')
    b1.set_ylabel('i') 
    b1.set_title('All Up < Tc')
    
    b3 = plt.subplot2grid((2,2),(0,1))
    b3.imshow(cf2, extent=extent, origin='upper') 
    b3.set_xlabel('j')
    b3.set_ylabel('i') 
    b3.set_title('All Up > Tc')

    b2 = plt.subplot2grid((2,2),(1,0))
    b2.imshow(cf3, extent=extent, origin='upper') 
    b2.set_xlabel('j')
    b2.set_ylabel('i')  
    b2.set_title('Random < Tc')

    b4 = plt.subplot2grid((2,2),(1,1))
    b4.imshow(cf4, extent=extent, origin='upper') 
    b4.set_xlabel('j')
    b4.set_ylabel('i')  
    b4.set_title('Random > Tc')
    plt.tight_layout()
    plt.show()

###### Auto Corelation ##########
def autoc(mat):       
    fftx = np.fft.fft2(mat)
    ret = np.fft.ifft2(fftx * np.conjugate(fftx))
    ret = np.fft.fftshift(ret)
    ret = np.real(ret)
    return ret

###############################################
#Part 13

def Mvsh(start, stop, J, Beta, Tc, Ti, Tf, hi1, hf, hstep, n, po):
    its1 = []
    M1 = []
    Ti1 = Ti
    h1 = []
    #loop stuff
    hits = 0
    mi = 0
#make matrix
    mup = np.ones((n,n))  #create nxnarray of 1's
    mdwn = np.negative(np.ones((n,n)))  #create nxnarray of 1's
    mrdm = np.sign(np.random.random((n,n))-po)

#run MC for up matrix
    while True: # matrix loop
        Ti = Ti1
        zT=0
        if mi == 1:
            mat1 = np.copy(mdwn)
            Title = "All Down"
        elif mi == 2:
            mat1 = np.copy(mrdm) 
            Title = "Random"
        else:
            mat1 = np.copy(mup)  
            Title = "All Up"
            
        while True: # Temp loop
            hi = hi1
            M1=[]
            h1 =[]
            while True: # hloop
                B = Beta / (Ti)
                itsa, Ma, m = MCflip_Mag(J,hi,mat1,B,start,stop) 
                M1.append(Ma)
                #its1.append(itsa)
                #mats.append(mat1)
                hi += hstep
                h1.append(hi)
                if hi >= hf:
                    break
            if Ti == Tf:
                dev2 = np.std(M1,axis=1).tolist()       #creats array of stddevs
                #print dev2
                Mean2 = np.mean(M1,axis =1).tolist()    #creates array of M means
                #print Mean2
                h2=h1
                #print "h2 %s" %h2

            else:
                dev1 = np.std(M1,axis=1).tolist()       #creats array of stddevs
                #print dev1
                Mean1 = np.mean(M1,axis =1).tolist()    #creates array of M means
                #print Mean1
                h=h1
                #print h
                
            Ti = Tf
            zT += 1
            if zT > 1:
                break
        print mi
        ax = plt.subplot2grid((3,1),(mi,0))
        ax.errorbar(h, Mean1, yerr=dev1, label='T<Tc' ) 
        ax.errorbar(h2, Mean2, yerr=dev2, label='T>Tc' )             
        ax.set_xlabel('h')
        ax.set_ylabel('<M>')
        ax.set_title(Title)
        ax.legend(loc=1) 

        mi+=1
        if mi == 3:
            break
    plt.show()

#show diffrent ways to create a random matrix
def creatematrix():
    mat = np.sign(np.random.random((n,n))-po) # define and initialize the matrix or could use rand(n,n)
    mat2 = np.sign(np.random.random(n**2)-po).reshape(n,n)   # or use reshape

#Set Variables
start = 100
stop = 50000
J = 1.0                                 #Exchange Energy [Units??]
h = 0.0                                 #External applied magnetic field [Units??]
B = np.log(1+np.sqrt(2))/(2*J)      #Beta [Units?]
po = 0.5    #probability of spin up
Tc = (2*J)/np.log(1+np.sqrt(2))

#P1_13
Mvsh(5000, 10000, J, B,Tc, 0.1, 2, 0.05, 10, 0.05, 50, po)       #start, stop, J, h, Beta, Tc, Ti(%Tc), Tf(%Tc), 
                                                                        #hi, hf, hstep, n, po

#P1_10-12
#MvsT(5000, 10000, J, h, B, Tc, 0.05, 3, 0.05, 50, po)      #start, stop, J, h, B ,Tc, Ti(%Tc), Tf(%Tc), Tstep, 
                                                    #n(Matrix size), po(probability fo randmat) -rand/up matrix made in func
#P1_9
#mat =np.ones((50,50))                         #create nxn array of 1's  
#run the desired MC and plots it
#test_con(100, 20000, J, h, B, mat)     #start, step, stop, J, h, B, mat #converges around 10000

#P1_8
#mat = np.array([[1,1,1],[1,1,1],[1,1,-1]]) #for delta E Test
#test_delta_E(J,h,mat)                      #it works!!!

#print wrap_i(0,2,-1)
#print wrap_i(1,2,-1)
#print wrap_i(0,2,1)
#print wrap_i(1,2,1)


