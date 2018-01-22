# -*- coding: utf-8 -*-
import math
import numpy as np
import matplotlib.pyplot as plt

i = 0
C = 0 #Counts in Circle
PI=[]
its = []
it=100
itmax = 100000    #iterations
while True:
    while True:
        x = np.random.random() 
        y = np.random.random()
        R=math.sqrt(x**2 + y**2)
        if R <= 1:
            C += 1.0
        i += 1
        if i == it:
            break
    PI.append(C/i*4)
    its.append(it)
    it +=100
    if it >= itmax:
        break

plt.plot(its, PI, 'r.-')
plt.xlabel('Iterations')
plt.ylabel(r'Value of $\pi$')  #π

plt.show()
#fig=plt.figure(figsize=(6,6), dpi=200) 

#ax1=fig.add_subplot(111)
#plt.grid(True)

