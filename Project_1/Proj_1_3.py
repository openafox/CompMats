
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
import time

#KMC 
def KMC(n, step_time, stop_time, its): #n= mat size, time(s) #for future nv = number of vacancies
    nv = 1      #1 vancancy
    v = np.random.randint(n**2-1)         #randomly pick vacancy position
    print "v: %d" %v
    Ind = np.arange(n**2)             #make array 0 - n^2
    Ad_Bk = []
    disp_sq_v = []
    disp_sq_a = []

    it = 0 
    i = 0
    #set up ploting
    #mx = matplt_setup(n)

    ########################
    #create address book
    while True:
        Ad_Bk.append(list(np.unravel_index(Ind[i],(n,n))))          #create adress book np.unravel converts index to i,j 
        i += 1                                                      #ex. 5 in a 3x3 is (1,2)
        if i == n**2:
            break
    mat = Ind.reshape(n,n)      #make matrix of numbers
    #make copies to maintain
    mat_i = mat                     
    Ad_Bk_i = Ad_Bk
    ###############################
    #itteration loop
    while True:              
        mat = mat_i                #reset matrix and adbook
        Ad_Bk = Ad_Bk_i
        disp_sq = np.array([])
        r=[]
        i = 0
        Time = 0
        T = []
        while True:             #define rate vector
            r.append(0.25)
            i += 1
            if i == nv*4:
                break
        rt = np.sum(r)          #get rate total
        #print "r: %s" %r
        ########################
        #simulation loop
        while True:         
            r_rand = np.random.random()    #pick random rate (rand 0-rt)
            #t = -np.log(np.random.random())/rt  #select time step
            Time += 1
            #print "Time: %d" %Time

            x = rt*r_rand           #define x
            #print "x: %s" %x
            #rate selector
            i = 0
            while True:     
                if r[i] >= x:
                    break
                else:
                    x = x-r[i]
                    i += 1
            ###########################
            #move based on selected rate            
            #calling r[i%4=0]=North(i-1), r[i%4=1]=South(i+1), r[i%4=2]=West(j-1), r[i%4=3]=East(j+1)
            norm = i        #(i - 4) * i//4     #normalize to 0-4 ie which vacancy (at this point does not matter but will need a vacancy list...
            #print norm
            #determin (-) or (+) shift i.e k = +-1
            k = (((norm+1)*2)-5)       
            k = k/np.abs(k)  

            #get initial vacancy positions
            i_i = wrap_i(Ad_Bk[v][0],n)
            j_i = wrap_i(Ad_Bk[v][1],n)
            if norm%2 == 0:
                #print "i %s" %k
                #get final vacancy positions
                i_f = wrap_i(Ad_Bk[v][0] + k, n)
                j_f = wrap_i(Ad_Bk[v][1], n)
                g = mat[i_f][j_f]
                #print Ad_Bk[v], i_i, j_i, i_f, j_f
                #edit adress book
                Ad_Bk[v][0] += k
                Ad_Bk[g][0] += -k

            else:
                #print "j %s" %k
                #get final vacancy positions
                i_f = wrap_i(Ad_Bk[v][0], n)
                j_f = wrap_i(Ad_Bk[v][1] + k, n)
                g = mat[i_f][j_f]
                #print Ad_Bk[v], i_i, j_i, i_f, j_f
                #edit adress book
                Ad_Bk[v][1] += k
                Ad_Bk[g][1] += -k
            #edit matrix
            mat[i_f][j_f] = mat[i_i][j_i]
            mat[i_i][j_i] = g
            ####################################

            #####change to np.arrays for manipulation
            Ad_Bk = np.array(Ad_Bk)
            Ad_Bk_i = np.array(Ad_Bk_i)

            if Time % step_time ==0 :
            #calc displacement matrix
                #print Ad_Bk_i
                #print "Ad_Bk_i[:,0]: %s" % Ad_Bk_i[:,0]
                i_dist = np.subtract(Ad_Bk_i[:,0], Ad_Bk[:,0])
                #print "i_dist: %s" %i_dist
                j_dist = np.subtract(Ad_Bk_i[:,1], Ad_Bk[:,1])
                if len(T)>0:
                    disp_sq = np.vstack((disp_sq,np.add(np.square(i_dist),np.square(j_dist))))
                else:
                    disp_sq = np.add(np.square(i_dist),np.square(j_dist))
                #print "disp_sq: %s" %disp_sq
                
                T.append(Time)

            #Plot live matrix
            #update_mat(mx,mat, v)       #for error checking

            if Time >= stop_time:
                break
        ##########################
        disp_sq_d = np.copy(disp_sq)
        #print "disp_sq: %s" %disp_sq

        disp_sq_d = np.delete(disp_sq_d,v,axis=1).tolist()
        #print "disp_sq_d: %s" %disp_sq_d
        if it == 0:
            disp_sq_v = disp_sq[:,v]
            disp_sq_a = np.sum(disp_sq_d, axis=1)


        else:

            disp_sq_v = np.vstack((disp_sq_v,disp_sq[:,v]))
            disp_sq_a = np.vstack((disp_sq_a, np.sum(disp_sq_d, axis=1)))
        #print "disp_sq_v_%s: %s" %(it, disp_sq_v)
        #print "disp_sq_a_%s: %s" %(it, disp_sq_a)
        
        it += 1
        if it >= its:
            break
    ##########################

    T=np.array(T)
    #print "T: %s" %T
    dev_v = np.std(disp_sq_v,axis=0)      #creats array of stddevs
    #print "dev_v: %s" %dev_v
    Mean_v = np.mean(disp_sq_v,axis =0)    #creates array of M means
    #print "Mean_v: %s" %Mean_v
    dev_a = np.std(disp_sq_a,axis=0)      #creats array of stddevs
    #print "dev_a: %s" %dev_a
    Mean_a = np.mean(disp_sq_a,axis =0)#/(n**2)    #creates array of M means
    #print "Mean_a: %s" %Mean_a

    plt.errorbar(T, Mean_v,yerr=dev_v, label="Vacancy") 
    plt.errorbar(T, Mean_a,yerr=dev_a, label="Atoms") 
    plt.xlabel('t')
    plt.ylabel('<d2>')
    plt.legend(loc=2)
    plt.show()



def update_mat(mx, mat, v):
    mat_up = np.copy(mat)
    mat_up[mat_up != v] = 1
    mat_up[mat_up == v] = 0
    mx.set_data(mat_up)
    plt.draw()
    #time delay to see motion
    time.sleep(0.2)  

def update_plt(hl, new_data_x, new_data_y):
    hl.set_xdata(numpy.append(hl.get_xdata(), new_data_x))
    hl.set_ydata(numpy.append(hl.get_ydata(), new_data_y))
    hl.relim()
    hl.autoscale_view()
    plt.draw()
def matplt_setup(n):
    #Plotting   (Matricies can be plotted with plt.matshow(mat) but I like it to be cleaner so extra stuff
    extent = [-0.5, n-0.5, n-0.5, -0.5]       #use shape to calc the "extent"

    mx = plt.imshow([[]], extent=extent, origin='upper', cmap=plt.get_cmap('hot'), interpolation='nearest', vmin=0, vmax=1)    #  plt.subplot2grid((1,1),(0,0))          #Plot it in the upper left coner
    plt.ion()
    plt.show()

#plt.axis([0, 1000, 0, 1])

#for i in range(1000):
   # y = np.random.random()
    #plt.scatter(i, y)
    #plt.draw()
    #time.sleep(0.05)


    #mx.xlabel('j')
    #mx.ylabel('i')
    #hl, = plt.plot([], [])          #set up plot
    #mx = plt.matshow(np.ones((n,n)))
    return mx



#Function to Wrap an Index in either direction
def wrap_i(i,n):            #(-)direction works automaticly in python     
    i = (i)%(n)
    return i



KMC(8, 100, 1000, 1000)    #(n, step_time, stop_time, its):
