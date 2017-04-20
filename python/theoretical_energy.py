# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 19:41:45 2017

@author: johan
"""
import numpy as np

def half_energy(tau,d,beta,hw,part_type):
    if(tau==0):
        hwb = hw*beta
        frac1 = (2*np.sinh(hwb/2)**2/np.sinh(hwb))**d
        frac2 = np.tanh(hwb/2)**(-1)
        frac3 = np.tanh(hwb)**(-1) 
    else:
        P = np.round(beta/tau)
        hwb = hw*beta
        x = hwb/P
        b = 1 + 0.5*x**2 + 0.5*x*np.sqrt(4+x**2)
        bder = x + 0.5*np.sqrt(4+x**2) + 0.5*x**2/np.sqrt(4+x**2)
        bfrac = bder/b
        frac1 = ((b**P-1)**2/(b**(2*P)-1))**d
        frac2 = (b**P+1)/(b**P-1)
        frac3 = (b**(2*P)+1)/(b**(2*P)-1)
    
    E=0
    if(part_type=='dis'):
        E = 2*0.5*hw*d*(0.5+1/(np.exp(hwb)-1))
    elif(part_type=='bos'):
        E = 0.5*hw*d*(frac2 + frac3*frac1)/(1+frac1)
    elif(part_type=='fer'):
        E = 0.5*hw*d*(frac2 - frac3*frac1)/(1-frac1)
    if(tau!=0):
        E *= bfrac
    return E