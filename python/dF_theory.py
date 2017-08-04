# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:22:24 2017

@author: johan
"""

import numpy as np
import matplotlib.pyplot as plt

hw = 3.0
beta=1.0
m = 1.0
n = 50
m2 = np.linspace(1.0,10.0,n)


exp1 = (np.exp(2*hw*beta*np.sqrt(m/m2))-1)/(np.exp(2*hw*beta)-1)
exp2 = ((np.exp(hw*beta*np.sqrt(m/m2))-1)/(np.exp(hw*beta)-1))**2

dF1 = -np.log(exp1)/beta
dF2 = -np.log(exp2)/beta

expm1 = (np.exp(hw*beta)-1)**2/(np.exp(hw*beta*2)-1)
expm2 = (np.exp(hw*beta*np.sqrt(m/m2))-1)**2/(np.exp(hw*beta*2*np.sqrt(m/m2))-1)

e1 = (1-expm1)/(1+expm1)
e2 = (1-expm2)/(1+expm2)
dFm1 = -np.log(expm1)/beta*np.ones(n)
dFm2 = -np.log(expm2)/beta


plt.figure(1)
plt.clf()
#plt.plot(m2,dF2)
plt.plot(m2,dF1-dF2,label='mass contribution')
plt.plot(m2,dFm1,label='dF m1')
plt.plot(m2,dFm2,label='dF m2')
plt.plot(m2,dFm2+dF1-dF2,label='dF m2 + mass')
plt.legend(loc='lower right',fontsize=20)
#print(dF1)
#print(dF2)
#print(dF1-dF2)


T = np.linspace(1,50,100)
beta = 11.6/T
ZP = np.exp(hw*beta/2)/(np.exp(hw*beta)-1)
err = np.log(1-np.exp(-hw*beta)+np.exp(-2*hw*beta))/beta

plt.figure(2)
plt.clf()
plt.plot(T,ZP)
plt.plot(T,-np.log(ZP)/beta)
plt.plot(T,err)