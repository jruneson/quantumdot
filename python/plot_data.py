import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams



datasets = [['Distinguishable','run4'],['Bosons','run5'],['Fermions','run6']]
settings = datasets[2]

folder = settings[1]
tasks = [('dist','../'+folder+'/results.dat')]

plt.figure(1)
plt.clf()
plt.rc('text',usetex=True)
#plt.rc('font',family='serif')
rcParams.update({'font.size': 22})

d=2
hw = 3.0
tau = 0.15

kB = 1/11.6045


for t in tasks:
    name = t[0]
    data = np.loadtxt(t[1],comments='%')
    beta = data[:,1]
    T = 1./(kB*beta)
    epot = data[:,2]/hw
    epot_err = data[:,3]/hw
    plt.figure(1)
    #plt.plot(T,epot,'x-',label=name)   
    plt.errorbar(T,epot,epot_err,marker='x',label='Potential energy')
    
    ekin = data[:,4]/hw
    ekin_err = data[:,5]/hw
    plt.errorbar(T,ekin,ekin_err,marker='x',label='Kinetic energy')
    
    evir = data[:,6]/hw
    evir_err = data[:,7]/hw
    plt.errorbar(T,evir,evir_err,marker='x',label='Virial kinetic')

    T = np.linspace(0.1,80,100)
    beta = 1.0/(kB*T)
    P = np.round(beta/tau)
    hwb = hw*beta
    x = hwb/P
    b = 1 + 0.5*x**2 + 0.*x*np.sqrt(4+x**2)
    bder = x + 0.5*np.sqrt(4+x**2) + 0.5*x**2/np.sqrt(4+x**2)
    frac1 = (2*np.sinh(hwb/2)**2/np.sinh(hwb))**d
    frac2 = np.tanh(hwb/2)**(-1)
    frac3 = np.tanh(hwb)**(-1)
    E = 2*0.5*d*(0.5+1/(np.exp(hwb)-1))
    Eb = 0.5*d*(frac2 + frac3*frac1)/(1+frac1)
    Ef = 0.5*d*(frac2 - frac3*frac1)/(1-frac1)
    plt.plot(T,E,'k--')
    plt.plot(T,Eb,'k--')
    plt.plot(T,Ef,'k--')
    
    plt.ylim([0,5])
    plt.xlabel('$T~(\mathrm{K})$')
    plt.ylabel('$\mathrm{Energy}~(\hbar\omega_0)$')
    plt.legend(loc='lower right',fontsize=20)
    plt.title(settings[0])
    
plt.show()

