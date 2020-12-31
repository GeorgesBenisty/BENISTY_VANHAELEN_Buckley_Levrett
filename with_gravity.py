# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 12:09:50 2020

@author: Benisty Georges, Luke Vanhaelen


In this file, all codes for the petroleum recovery with gravity
are put here in one file !
We have :
    - Upwind
    - Lax-Wendroff
    - Godunov

Common parameters are initialized at the beginning
"""

# Variables :
import numpy as np
import matplotlib.pyplot as plt
a,b=0,2000 #Beginning and end of the space model
alpha=0.1472 #Initial conditions
N=200 #points
h=(b-a)/N
X=np.linspace(a,b,N) #Generate N points between a and b
T=3 #boucles de temps (secondes ?)
k=h/6 #lié à la CFL
r=k/h # dx/dt
S1=np.zeros(len(X)) #init saturation vector for upwind
S2=S1 # for LW
S3=S1 # for Godunov
S1[0]=1 #boundary conditions
S2[0]=S1[0]
S3[0]=S1[0]
U=np.zeros(N) #for Godunov
U[0]=1 #Initial conditions also apply !

K=0.1 #permeabilite
g=9.81 #gravity acceleration, m/s²
S1=np.zeros(len(X)) #init saturation vector
S1[0]=1 #boundary conditions
rho_w,rho_o = 0.2070,0.1520 #density
phi=0.4 #porosity
beta = (rho_w - rho_o)*(K/phi)*g #correction factor due to gravity

#Functions :
def fwater(s):
    return s**2

def foil(s):
    return (1-s)**2/4


def G(a,b):
    if -alpha+(beta*fwater(a))<=0:
        G=(fwater(a)*(alpha+beta*foil(a)))/(fwater(a)+foil(a))
    else:
        G=(fwater(a)*(alpha+beta*foil(b)))/(fwater(a)+foil(b))
    return(G)

#Upwind :
for n in range(0,T-1):
    for i in range(1,len(S1)-1):
       a=S1[i]
       b=(k/h)
       c=G(S1[i],S1[i+1])-G(S1[i-1],S1[i])
       S1[i] = a - (b*c)
    
#Lax-Wendroff :
for n in range(0,T-1):
    for i in range(1,len(S2)-1):
       S2[i] = (0.5*(S2[i+1]+S1[i-1]))-(0.5*r*(G(S2[i+1],S2[i])-G(S2[i-1],S2[i])))

#Godunov :
for n in range(0,T-1):
    for i in range(1,len(S3)-1):
        d=S3[i]-k*(G(S3[i],S3[i+1])-G(S3[i-1],S3[i])) #=> U_m^(1)
        e=S3[i-1]-k*(G(S3[i-1],S3[i])-G(S3[i-2],S3[i-1])) #=> U_m-1^(1)
        a=S3[i+1]-k*(G(S3[i+1],S3[i])-G(S3[i],S3[i+1])) #=> f(U_m-1^(1))
        S3[i] = 0.5*(S3[i]+d-k*(G(d,a)-G(e,d)))

plt.clf()
plt.plot(X,S1,label='Upwind scheme')
plt.plot(X,S2,label='Lax-Friederich scheme')
plt.plot(X,S3,label='Godunov scheme')
plt.xlabel('length (m) on z-direction')
plt.ylabel('saturation of the shale layer')
plt.title(
'Petroleum recovery with gravity \n Time = '+str(T)+' s, alpha = '+str(alpha))
plt.legend(loc='best')
