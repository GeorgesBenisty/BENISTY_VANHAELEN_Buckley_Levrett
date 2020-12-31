# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 18:31:11 2020

@author: Benisty Georges
"""

##########################" TP4 Ma411 #################"

"""
In this file, all codes for the petroleum recovery without gravity
are put here in one file !
We have :
    - Upwind
    - Lax-Wendroff
    - Godunov

Common parameters are initialized at the beginning
"""

# Common parameters :
#
# Variables :
import numpy as np
import matplotlib.pyplot as plt
a,b=0,2000 #Beginning and end of the space model
alpha=0.1472 #Initial conditions
N=200 #points
h=(b-a)/N
X=np.linspace(a,b,N) #Generate N points between a and b
T=180 #boucles de temps (secondes ?)
k=h/3 #lié à la CFL
r=k/h # dx/dt
S1=np.zeros((len(X),T)) #init saturation vector for upwind
S2=S1 # for LW
S3=S1 # for Godunov
S1[0,:]=1 #boundary conditions
S2[0,:]=S1[0,:]
S3[0,:]=S1[0,:]
U=np.zeros(N) #for Godunov
U[0]=1 #Initial conditions also apply !

#Functions :
def fwater(s):
    return s**2

def foil(s):
    return (1-s)**2/4

def f(s):
    return alpha*fwater(s)/(fwater(s)+foil(s))


# Loop to store values of saturation :

#Upwind :
for n in range(0,T-1):
    for i in range(1,len(S1)-1):
       S1[i,n+1] = S1[i,n] - r*(f(S1[i,n])-f(S1[i-1,n]))

#Lax-Friederich:
for n in range(0,T-1):
    for i in range(1,len(S2)-1):
       S2[i,n+1] = (0.5*(S2[i+1,n]+S2[i-1,n]))-(0.5*r*(f(S2[i+1,n])-f(S2[i-1,n])))

#Godunov:
for n in range(0,T-1):
    for i in range(1,len(S3)-1):
        U[i]=S3[i,n]-r*(f(S3[i,n])-f(S3[i-1,n]))
        S3[i,n+1] = 0.5*(S3[i,n]+U[i]-r*(f(U[i])-f(U[i-1])))


# Plotting the curve :
plt.clf()
plt.plot(X,S1[:,-1],label='Upwind scheme')
plt.plot(X,S2[:,-1],label='Lax-Friederich scheme')
plt.plot(X,S3[:,-1],label='Godunov scheme')
plt.xlabel('length (m) on z-direction')
plt.ylabel('saturation of the shale layer')
plt.title(
'Petroleum recovery without gravity \n Time='+str(T)+' s, '+'alpha ='+str(alpha)
)
plt.legend(loc='best')