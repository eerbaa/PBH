# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 11:57:54 2023
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import scipy.interpolate as spip

#Parameters

m = 0.01
alpha = 1

#Framework

def mesh_r(start, stop, num):
    return(np.linspace(start,stop, num))

#Functions

def Th(t):
    return(100) #to change later

def Temp(T0, mesh):
    return((spip.interp1d(mesh, T0)))

def g(Tem):
    return(1) #to change later 

def E(mesh, t, T0):
    T = Temp(T0, mesh)
    def dEdr(r, E):
        return((-alpha**2)*np.sqrt(((T(r)**3)*((E+m)**3))/((E*2)+m*E))*(1-(m/(E+m))**2)**(2))
    y0 = np.array([Th(t)])
    sol = spi.solve_ivp(dEdr,[0,mesh[mesh.size-1]], y0, 'RK45', mesh)
    res = sol.y
    res = np.reshape(res, res.size)
    return(res)

   
#test
num =20 #number of points in the mesh
mesh = mesh_r(0, 180000000,num)
t = 0
T03 = 0.0001*np.ones(num)
T01 = np.ones((1,num))
T02 = 0.1*np.ones((1,num))
#T = Temp(T0, mesh)
#plt.plot(mesh, T(mesh))
#print(T(2.6))

#plot

E1 = E(mesh,t, T03)
E2 = E(mesh, t, T01)
E3 = E(mesh, t, T02)
plt.xlabel('r')
plt.ylabel('E')
plt.plot(mesh,E1, label = 'T0 = 0.0001')
#plt.plot(mesh,np.transpose(E2), label = 'T0 = 1')
#plt.plot(mesh,np.transpose(E3), label = 'T0 = 0.1')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
