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
#hw = 3.0/7300
hw = 3.0
hbar = 0.6582
beta = 1.0
T = 1.0/(kB*beta)
E = hw/10.0
#m = 4.002602*1.6605e-27
m_rel = 1000
m = 9.1e-31*m_rel

#eps = 0.8804
eps = 1.0
#s = 0.2556e-9/5.29177e-11
s = 50.0
r = np.linspace(0.8*s,2*s,1000)


def V(r):
    return 4*eps*((s/r)**12-(s/r)**6)



"""plt.figure(1)
plt.clf()
plt.plot(r,V(r))
plt.plot([r.min(),r.max()],[E, E],'k--')
plt.xlabel('$r$')
plt.ylabel('$V(r)$')"""

tol = 0.0001
b = s+1.0
a = s/10.0
rmin = s
while(abs(V(rmin)-E)>tol or V(rmin)-E>0):
    rmin = (a+b)/2.0
    diff = V(rmin)-E
    if(diff<0):
        b = rmin
    elif(diff>0):
        a = rmin

rmax = s+1.0
r = np.linspace(rmin, rmax, 10000)
integrand = 1.0/np.sqrt(2*(E-V(r))*5.72e-26/m)
"""plt.figure(2)
plt.clf()
plt.plot(r,integrand)"""
dr = r[1]-r[0]
print(sum(integrand)*dr/20)
print(hbar/hw *2*np.pi /50)
print(np.sqrt(1.0/(3.67e-5*hw*m_rel)))
