# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 14:27:36 2017

@author: johan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.rc('text',usetex=True)
plt.rc('font',family='helvetica')
rcParams.update({'font.size': 22})

kB = 1.0/11.6045
m_rel = 0.067
hw = 1.2 #/m_rel
hbar = 0.6582
beta = 1.0#*m_rel
P=50
T = 1.30/(kB*beta)
E = hw*3*2/P

#m = 4.002602*1.6605e-27
m = 9.1e-31*m_rel

#eps = 0.8804
#eps = hw*0.1
#s = 0.2556e-9/5.29177e-11
l0 = 1500
r = np.linspace(0.1*l0,2*l0,1000)


def V(r):
    return 2109.47/(r*P)
    #return 4*eps*((s/r)**12-(s/r)**6)/P



plt.figure(1)
plt.clf()
plt.plot(r,V(r))
plt.plot([r.min(),r.max()],[E, E],'k--')
plt.xlabel('$r$')
plt.ylabel('$V(r)$')

tol = 0.0001
b = l0
a = l0/20.0
rmin = l0
while(abs(V(rmin)-E)>tol or V(rmin)-E>0):
    rmin = (a+b)/2.0
    diff = V(rmin)-E
    if(diff<0):
        b = rmin
    elif(diff>=0):
        a = rmin

print(rmin)
rmax = rmin+2
r = np.linspace(rmin, rmax, 10000)
integrand = 1.0/np.sqrt(2*(E-V(r))*5.72e-26/m)
plt.figure(2)
plt.clf()
plt.plot(r,integrand)
dr = r[1]-r[0]
dt_LJ=sum(integrand)*dr/20
dt_other=hbar/hw *2*np.pi /np.sqrt(1+4*P/(hw*hw*beta*beta))/20
print(dt_LJ)
print(dt_other)
print(dt_other/dt_LJ)
#print(np.sqrt(1.0/(3.67e-5*hw*m_rel)))
