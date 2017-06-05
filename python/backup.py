# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:47:16 2017

@author: johan
"""

def plot_energies_vs_t(f,fig_nr,n=200000,P=20,conn=False):
    #data = np.genfromtxt(f+'Pot_energy.dat',max_rows=n)
    """epot = data[:,1]
    data = np.genfromtxt(f+'Kinetic_energy.dat',max_rows=n)
    ekin = data[:,1]    
    data = np.genfromtxt(f+'Kin_en_virial.dat',max_rows=n)
    evir = data[:,1]"""
    #data = load_lines(f+'Total_energy_cl.dat',n)
    data = np.genfromtxt(f+'Total_energy.dat',max_rows=n)
    #data = np.genfromtxt(f+'X_coord_n1p1.dat',max_rows=n)
    signal = data[:,1]
    #data = np.genfromtxt(f+'Total_energy_cl.dat',max_rows=n)
    #etotc= data[:,1]
    t = data[:,0]
    """data = np.genfromtxt(f+'Pot_energy_cl.dat',max_rows=n)
    epotc = data[:,1]
    data = np.genfromtxt(f+'Kin_energy_cl.dat',max_rows=n)
    ekinc = data[:,1]
    data = np.genfromtxt(f+'Twopart_energy.dat',max_rows=n)
    eintc = data[:,1]
    data = np.genfromtxt(f+'Spring_pot_en.dat',max_rows=n)
    esprc = data[:,1]"""
    """plt.figure(fig_nr)
    plt.clf()
    plt.plot(t,etot,marker='.',label='etot')"""
    plt.figure(fig_nr)
    plt.clf()
    dt = t[1]-t[0]
    N = len(signal)
    xf = np.linspace(0.0,1.0/(2.0*dt),N/2)
    yf = np.fft.fft(signal)*2.0/N
    print(N)
    plt.plot(xf,np.abs(yf[:N//2]))
    plt.xlim([0,15])
    plt.ylim([0,6])
    plt.xlabel(r'$f~(\mathrm{ps}^{-1})$')
    plt.ylabel(r'$\mathrm{FFT~of~} E$')
    plt.ylabel(r'$\mathrm{FFT~of~}\langle WE\rangle /\langle W \rangle$')
    """w0 = 3/0.6582/np.sqrt(P)
    k = (2*20/(3*2))**2
    P2 = P*(1+conn)
    j = np.linspace(1,P2,P2)
    ws = w0*np.sqrt(1+k*np.sin(np.pi*j/P2)**2)"""
    #plt.plot(ws/(2*np.pi),0.8*np.ones(len(ws)),'x')
    #plt.ylim([0,0.5])
    #plt.plot(t,epot,marker='.',label='epot')
    #plt.plot(t,ekin,marker='.',label='ekin')
    #plt.plot(t,evir,marker='.',label='evir')
    #plt.plot(t,etotc,marker='.',label='etotc')
    #plt.plot(t,esprc,marker='.')
    #plt.xlim([0,10])
    #plt.plot(t,ekin)
    #plt.plot(t,evir)
    
    
"""    
def plot_rAB_dist(settings,ax=None):
    f = settings[0]
    d = settings[1]
    name = settings[2]
    marker = settings[4]
    color= settings[5]
    #f = f_name[0]
    #name = f_name[1]
    data = np.loadtxt(f+'Pair_correlation.dat')
    r = data[:,0]
    p = data[:,1]
    p_err = data[:,2]
    p1 = data[:,3]
    p1_err = data[:,4]
    p2 = data[:,5]
    p2_err = data[:,6]
    a = np.sqrt(1.0/(3.67e-5*3))
    #print(a)
    #gaussian = np.exp(-(r/a)**2)
    p2part = np.exp(-0.5*(r/a)**2)
    #normalize(gaussian,0,r,d)
    normalize(p2part,0,r,1,1)
    normalize(p,p_err,r,d)
    p *= 1e6
    p_err *= 1e6
    #normalize(p1,p1_err,r,d)
    #normalize(p2,p2_err,r,d)
    #p /= np.pi * ((r+dr)**2-r**2)
    #p_err /= np.pi * ((r+dr)**2-r**2)
    #plt.plot(r,np.log(p))
    if ax is None:
        ax = plt.gca()
    ax.plot(r,p,color=color)
    ax.errorbar(r[::3],p[::3],p_err[::3],linestyle='None',marker=marker,color=color)
    #plt.errorbar(r,p1,p1_err,label='$p_1(r)$')
    #plt.errorbar(r,p2,p2_err,label='$p_2(r)$')
    #plt.plot(r,gaussian,label='theory $p_n(r)$')
    
    #plt.plot(r,p2part/170,label='theory $p(r_\mathrm{AB})$')
    #plt.legend(loc='upper right',fontsize=22)
    #plt.title('Two fermions, Lennard--Jones-potential')"""

"""    
def load_lines(filename,num_lines):
    with open(filename) as file:
        head = [next(file) for x in range(num_lines)]
    for line in head:
        myarray = np.fromstring(line)
    print(head)
    return head"""
    
"""
def plot_rABs(fs):
    #fig = plt.figure(4)
    #plt.clf()
    #fig.set_size_inches(13,5.5,forward=True)
    #ax1 = fig.add_subplot(121)
    fig, (ax1,ax2) = plt.subplots( nrows=1, ncols=2, 
                             sharey=True,figsize=(13,5.5))
    #plt.subplot(121)    
    ax1.set_ylim([-0.5,2.1])
    ax1.set_title(r'$\mathrm{Without~metadynamics}$')
    #ax2 = fig.add_subplot(122, sharey=ax1)
    #plt.subplot(122)
    ax2.set_title(r'$\mathrm{With~metadynamics}$')
    ax2.set_ylim([-0.5,2.1])
    for settings in fs:
        if(settings[3]=='no MetaD'):
            plot_rAB_dist(settings,ax1)
        if(settings[3]=='with MetaD'):
            plot_rAB_dist(settings,ax2)
        if(settings[3]=='both'):
            plot_rAB_dist(settings,ax1)
            plot_rAB_dist(settings,ax2)
    ax1.set_xlabel('$r_\mathrm{AB}/a_0$')
    ax2.set_xlabel('$r_\mathrm{AB}/a_0$')
    ax1.set_ylabel('$p(r_\mathrm{AB})/10^{-6}$')
    ax1.annotate('$(a)$', xy = get_axis_limits(ax1))
    ax2.annotate('$(b)$', xy = get_axis_limits(ax2))
    plt.subplots_adjust(wspace=0)
    tmp = tic.MaxNLocator(8)
    ax1.xaxis.set_major_locator(tmp)"""