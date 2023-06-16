#!/usr/bin/env python
# coding: utf-8

# In[28]:


# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 11:57:54 2023
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import scipy.interpolate as spip

#Parameters

m = 0.1
alpha = 1
phi = 1

#Framework

def mesh_r(start, stop, num):
    return(np.linspace(start,stop, num))

#Functions

def Th(t):
    return(100) #to change later

def T(r, T0):
    return(spip.interp1d(np.transpose(T0), r))
   #return(1) #to change later

def g(Temp):
    return(1) #to change later 

def E(mesh, t, T0):
    def dEdr(r, E):
        return((-alpha**2)*np.sqrt(((T(r,T0)**3)*((E+m)**3))/((E*2)+m*E))*(1-(m/(E+m))**2)**2)
    y0 = np.array([Th(t)])
    sol = spi.solve_ivp(dEdr,[0,mesh[mesh.size-1]], y0, 'RK45', mesh)
    return(sol.y)

   
#test

mesh = mesh_r(0, 50,101)
t = 0
T0 = np.ones((1,101))
print()
#plot

E1 = E(mesh,t, T0)
plt.plot( mesh_r(0,30,100),np.transpose(E1))
plt.show()


# In[ ]:


def Te(T0, dt, t, mesh): #T0 is a np.array of length num
    if T0.size != num :
        raise Exception('The size of T0 should be num.')
    if num%2 == 0 :
        raise Exception('num must be odd')
    else
        r = mesh
        Tdt = np.zeros(1,num)
        dTdr = np.zeros(1, num) #we impose periodic boundary conditions 
        dTdr2 = np.zeros(1,num)
        dTdr2[0] = (T0[1]-2*T0[0]+T[num-1])/((r[1]-r[0])**2)
        dTdr2[num-1] = (T0[0]-2*T0[num-1]+T0[num-2])/((r[num-1]-r[num-2])**2)
        dTdr[1:num-1] = (T0[2:num]-T0[0:num-2])/(r[2:num]-r[0:num-2])
        dTdr2[1:num-1] = (T0[2:num]-2*T0[1:num-1]+T0[0:num-2])/(((r[2:num]-2*r[1:num-1]+r[0:num-2])/2)**2)
        E0 = np.zeros(1,num)
        E0[0:((num+1)/2)] = np.flip(E(rmin, rmax, num,t,T0[num/2])[0:(num+1)/2])
        E0[(num+1)/2:num] = E(rmin, rmax, num,t,T0[num/2])[1:(num+1)/2]
        Tdt = ((45/(((np.pi*alpha)**2)*T0*g(T0)))*dTdr +
               (45*(r**2)/(((np.pi*alpha*T0)**2)*g(T0)))*(dTdr)**2 +
               (30*(r**2)/(((np.pi*alpha)**2)*T0*g(T0)))*dTdr2 + 
               ((15*phi*(alpha)**2)/(8*((np.pi*T0)**3)*g(T0)))*np.sqrt(((T0*(E0+m))**3)/(E0**2 + m*E0))*((1-(m/(E+m))**2)**2)
              )*dt + T0
        return(Tdt)
    
#test

Tdt1 = Te()
            
            

