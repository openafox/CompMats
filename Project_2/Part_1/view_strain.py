import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def qrot(v1,ax,a):
    if len(np.shape(v1)) == 1:
        v1 = np.reshape(v1,(1, len(v1)))
    #if len(np.shape(ax)) == 1:
        #ax = np.reshape(ax,(len(ax),1))

    #print "ax: %s" %ax
    #print "v1: %s" %v1
# quaternion rotation of vectors v1 by angle a about axis axl through origin
    if np.dot(ax,ax) != 0:
        #print "v1: %s" %v1

        axl = ax/np.sqrt(np.dot(ax,ax))
        r = np.shape(v1)
        #print "r:", r

        if len(r) > 1:
            d  = r[0]
        else:
            d = 1

        zer = np.zeros((d,1))
        #print zer
        #if d == 1:
        #    zer = np.zeros((d,))
        #print "d: %s" %d
        dv = np.ones((d,1))
        #print "dv: %s" %dv
        aa = np.cos(a/2)*dv
        #print "aa: %s" %aa
        bb = np.sin(a/2)*(dv*axl)
        #print "bb: %s" %bb
        
        q  = np.hstack((aa,bb))
        #print "q: %s" %q
        #print "q: %s" %q
        #print "zer: %s" %zer
        v2a = np.hstack((zer,v1))
        #print "v2a: %s" %v2a
        rr = qmult(v2a,qconj(q))
        #print "rr: %s" %rr
        v2 = qmult(q,rr)
        v2 = v2[:,1:4]
    else:
        v2 = v1
    return v2

def qconj (q):
    qc = np.multiply(q, [1.0, -1.0, -1.0, -1.0])
    #print "qc: %s" %qc
    return qc

# Multiplication of a quaternion #np.dot(q1[1:],q2[1:])],\    

def qmult(q1,q2):
    #print "q1: %s" %q1
    #print "q2: %s" %q2
    qa = np.reshape(q1[:,0]*q2[:,0],(1,len(q1[:,0])))
    qb = -np.einsum('...j,''...j',q1[:,1:4], q2[:,1:4])  ####
    #print "qb: %s" %qb
    l = len(q1)
    #print "len: %s" %l
    qc = np.cross([q1[i,1:4] for i in range(l)],[q2[i,1:4] for i in range(l)])
    #print "qd2: %s" %q2[:,0]
    #print "ones: %s" % np.ones((1,3))
    qd = np.reshape(q2[:,0],(l,1))*np.ones((1,3))*q1[:,1:4]
    #print "qd: %s" %qd
    qe1 = q1[:,0]
    #print "qe1: %s" %qe1
    qe2 = q2[:,1:4]
    #print "qe2: %s" %qe2    
    qe = (np.reshape(q1[:,0],(l,1))*np.ones((1,3)))*q2[:,1:4]
    #print "qe: %s" %qe
    qab = (qa+qb)

    #print "qab: %s" %(qa+qb)
    #print "Qde: %s" %(qc + qd + qe)
    q3 = np.hstack((qab.T, qc + qd + qe))
    #print "q3: %s" %q3

    return q3

def view_strain(eps):

###  Set Up Figure ####
    font = {'family' : 'Helvetica',
            'weight' : 'normal',
            'size'   : 24}
    fig = plt.figure()
    plot = fig.add_subplot(111, projection='3d')
    #f = figure(12); clf;
    #axis image
    #set(gca,'box','off','Projection','perspective','FontName','Helvetica','FontSize',18)
    #set(gca,'LineWidth',[1.0]) + axis([-1,1,-1,1,-1,1]*1.5);
    plot.set_xlabel('x', **font)
    plot.set_ylabel('y', **font)
    plot.set_zlabel('z', **font)



    org = -1*np.array(([1,1,1]))
    vec1 = 2*np.array(([[1,0,0],[0,1,0],[0,0,1]]))

# Define path around box
    m = org
    m = np.vstack((m,m[-1] + vec1[0]))
    m = np.vstack((m,m[-1] + vec1[1]))
    m = np.vstack((m,m[-1] - vec1[0]))
    m = np.vstack((m,m[-1] - vec1[1]))
    m = np.vstack((m,m[-1] + vec1[2]))
    m = np.vstack((m,m[-1] + vec1[0]))
    m = np.vstack((m,m[-1] - vec1[2]))
    m = np.vstack((m,m[-1] + vec1[2]))
    m = np.vstack((m,m[-1] + vec1[1]))
    m = np.vstack((m,m[-1] - vec1[2]))
    m = np.vstack((m,m[-1] + vec1[2]))
    m = np.vstack((m,m[-1] - vec1[0]))
    m = np.vstack((m,m[-1] - vec1[2]))
    m = np.vstack((m,m[-1] + vec1[2]))
    m = np.vstack((m,m[-1] - vec1[1]))

   #print m 
    plot.plot(m[:,0],m[:,1],m[:,2], color='b', linewidth=2.0)
    
# strain the box
    m = m + np.dot(m, eps) #np.multiply(m, eps)
    plot.plot(m[:,0],m[:,1],m[:,2], color='g', linewidth=2.0)

# Rotate coordinate system
# First rotate so that the "a" vector is parallel to (1,0,0) 
    vec2 = vec1 + vec1 * eps
    ax   = np.cross(vec1[1],vec2[1])
    #print "ax: %s" %ax
    d1   = np.sqrt(np.dot(vec1[1],vec1[1]))
    d2   = np.sqrt(np.dot(vec2[1],vec2[1]))
    A3   = np.sqrt(np.dot(ax,ax))
    ang  = np.arcsin(A3/(d1*d2))

    #print "vec2: %s" %vec2
    vec3 = qrot(vec2,ax,-ang)
    #print "vec3: %s" %vec3
    m    = qrot(m,ax,-ang)
    #print "m: %s" %m

    #vec3 = qrot(vec2,ax,-ang);
    #m    = qrot(m,ax,-ang);


# Now rotate about (1,0,0) so that b vector is in x-y plane 
    v1 = vec3[1]
    #print vec3[1]
    #print v1
    v2 = np.array([0,1,1]) * v1
    #print v2
    v2 = v2/np.sqrt(np.dot(v2,v2))
    #print v2
    ang = np.arccos(np.dot(v2,[0,1,0]))
    #print "ang: %s" %ang
    ax  = [1,0,0];

    vec4 = qrot(vec3,ax,-ang)
    #print "vec4: %s" %vec4
    m    = qrot(m,ax,-ang)
    #print "m: %s" %m

    
    plot.plot(m[:,0],m[:,1],m[:,2], color='r', linewidth=2.0)

    plt.show()
############################

#a = 0.5
#mat = np.array([[1,0,0],[0,1,0],[0,0,1]])
#eps = a * mat
#print "eps: %s" %eps
#view_strain(eps)







