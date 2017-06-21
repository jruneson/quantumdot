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
    
    
    if 0:
        plot_energies(f3a,1,1,'tau',beta=1,label=r'$\beta=1~\mathrm{meV}^{-1}$')
        plot_energies(f3b,1,0,'tau',beta=2,label=r'$\beta=2~\mathrm{meV}^{-1}$',marker='o')
        #plt.plot([0,10],[2,2],'k--',label=r'$\mathrm{Theoretical~value}$')
        plt.legend(loc='lower right',fontsize=26)
        plt.title('Three particles')

    if 0:
        plot_theoretical_energies(1,True,2,3.0,'k',r'$\hbar\omega_0=3\,\mathrm{meV}$')
        plot_theoretical_energies(1,False,2,15.0,'b',r'$\hbar\omega_0=15\,\mathrm{meV}$')    
        plot_theoretical_energies(1,False,2,30.0,'r',r'$\hbar\omega_0=30\,\mathrm{meV}$')    
        plt.legend(loc='upper left',fontsize=24)
        plt.xlabel(r'$T\,(\mathrm{K})$')
        plt.ylabel(r'$E/\hbar\omega_0$')

    if 0:
        plot_energies(f3a,1,1,'P',label='beta1',linestyle='',marker='o',plot_all=0)
        plot_energies(f3b,1,0,'P',label='beta2',linestyle='',marker='D',plot_all=0)
        #plot_energies(fi3,1,0,'beta',label='Distinguishable',linestyle='',marker='s',plot_all=0)
        #plot_theoretical_energies(1,d=1)
        #plot_energies(f15,1,0,'P',label='Distinguishable',linestyle='-',marker='v',plot_all=0)
        #plot_energies(f22,1,0,'P',label='Unconnected, dt=1fs')
        #plt.plot([0,25],[16.4,16.4],'k--',label='Target value')
        #plt.plot([0,22],[3,3],'k--',label='Target value')
        plt.xlim([0,80])    
        plt.ylim([0,7])
        plt.legend(loc='upper left',fontsize=22)
        plt.title('1D')
        #plt.title(r'$\mathrm{Elliptic~QD,~}\tau=0.067~\mathrm{meV}^{-1}$')
        #plt.title(r'No Coulomb, $\hbar\omega_0=3\,\mathrm{meV}$')
        if 0:        
            plt.figure(2)
            plt.clf()
            plt.plot(x,euconn-econn,'o-')
            plt.plot(x,6-euconn,'o-')
            plt.grid(True)
        
        
    #plt.figure(3)
    #plt.title(r'$\mathrm{Without~Metadynamics}$')
    #free_energy_diff(fi3,fi4,4,beta=1,P=10)
    #plt.figure(4)
    #plt.title(r'$\mathrm{With~Metadynamics}$')
    #plt.title(r'Elliptic QD, $\gamma_\mathrm{screen}=1, P=15$')
    #plt.title(r'No Coulomb, $\hbar\omega=3\,\mathrm{meV}$')
    #plot_cv(f1,7,200000,opt='bdE')
    #plot_fes(f14,4)
    #data = np.loadtxt(f11+'Total_energy.dat')
    #plt.plot(data[:,0],data[:,2])
    #plot_2d_dist(fi2,3)
    #plot_1d_dist(fi,1,1)
        
    if 0:
        fig_nr = -1
        d = 1
        plt.figure(2)
        fig = plt.gcf()
        ax = plt.gca()
        r,pb,pberr=plot_rAB(fi3,fig_nr,1,d,'b','o',name='Boson',linestyle='-',show_errors=5)
        r,_,_=plot_rAB(fif,fig_nr,0,d,'g','D',name='Fermion',linestyle='-',show_errors=5)
        #r,pd,pderr=plot_rAB(fid,fig_nr,0,d,'r','s',name='Distinguishable',linestyle='-',show_errors=5)    
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        #pdth = plot_rAB_th(fig_nr,r,d,'dis')
        plt.legend(loc='upper right',fontsize=20)
        plt.title(r'$\mathrm{1D,~with~Metadynamics}$')
        #plt.ylim([-0.02,0.1])
        plt.ylim([0,0.2])
        #inset = np.zeros([50,50])
        #inset[0:50,0:50] = 
        #extent = [-3,4,-4,3]
        #ax.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.71,0.43,0.24,0.24])            
            rins = r[1:22:5]
            pins = pb[1:22:5]
            axins.errorbar(rins,pb[1:22:5],pberr[1:22:5],marker='o',color='b',linestyle='')        
            #axins.errorbar(rins,pd[1:30:5],pderr[1:30:5],marker='s',color='r',linestyle='')
            axins.plot(r[1:22],pbth[1:22],'b')
            #axins.plot(r[1:30],pdth[1:30],'r')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            plt.xticks(np.array([0,1,2]))        
            plt.yticks(np.array([0.15,0.16]))
            #plt.xticks(visible=False)
        #plt.yticks(visible=False)
        #mark_inset(ax,axins,loc1=1,loc2=4,fc="none",ec="0.5")
        #fig.sca(ax)
    
    if 0:
        plt.figure(6)
        plt.clf()
        plot_cont_sint(f18,4,50000,1.0,'FES')
        cvs,Fs = plot_cont_sint(f19,4,50000,1.0,'FES')
        plt.title('Blue = oo, Green = O')
        cvs = cvs[130:200]
        #plt.plot(cvs,cvs)

    if 0:
        beta = 1.0
        n = int(len(cvs))
        interval1 = range(100,145)
        interval2 = range(160,250)
        Fsmod = -Fs/beta-cvs
        poly1,cov1 = np.polyfit(cvs[interval1],Fsmod[interval1],1,cov=True)
        poly2,cov2 = np.polyfit(cvs[interval2],Fsmod[interval2],1,cov=True)
        cvs_ext = np.linspace(-25,0,100)
        cvs_ext2= np.linspace(0,40,100)
        plt.figure(5)
        plt.clf()
        plt.plot(cvs,Fsmod,'b')
        Fs_ext1 = np.polyval(poly1,cvs_ext)
        Fs_ext2=np.polyval(poly2,cvs_ext2)
        plt.plot(cvs_ext,Fs_ext1,'r')
        plt.plot(cvs_ext2,Fs_ext2,'g')
        plt.plot(cvs[interval1],Fsmod[interval1],'xr')
        plt.plot(cvs[interval2],Fsmod[interval2],'xg')
        plt.xlabel('$s$')
        plt.ylabel(r'$\log[p(s)e^{-s}]$')
        plt.ylim([-6,-2])
        plt.xlim([-30,10])
        print('slope1: {} +/- {}'.format(poly1[0],np.sqrt(cov1[0][0])))
        print('slope2: {} +/- {}'.format(poly2[0],np.sqrt(cov2[0][0])))
    if 0:
        slope1 = np.array([-1.4727,-1.2749, -1.1468,-1.0935,-1.0606,-1])
        s1err = np.array([0.0024,0.0017,0.0014,0.0015,0.0013,0.015])        
        slope2 = np.array([0.9540, 0.6830, 0.4774,0.3662,0.3007,0.1768])
        plt.figure(6)
        plt.clf()
        plt.errorbar(tau,-1-slope1,s1err,marker='.',linestyle='None',label='Slopes')
        poly3,cov3 = np.polyfit(tau,-1-slope1,2,cov=True)
        x = np.linspace(0,0.16)
        plt.plot(x,np.polyval(poly3,x),label='Quadratic fit')
        print('slope3: {} +/- {}'.format(poly3[0],np.sqrt(cov3[0][0])))
        plt.xlabel(r'$\tau~\mathrm{(meV)}^{-1}$')
        plt.ylabel(r'$\mathrm{-1-Slope1}$')
        plt.legend(loc='lower right', fontsize=20)
        polyexp,covexp = np.polyfit(tau,np.log(-1-slope1+1),1,cov=True)
        #plt.plot(x,np.exp(np.polyval(polyexp,x))/np.exp(polyexp[1])-1)
        #plt.plot(tau,-1-slope2,'o')
        #plot_cont_sint(f4,4,200000,'FES')
        #plt.title('bos')
    if 0:
        plot_cont(f4,5,1,200000,'bdE','etot')
        plt.figure(5)
        plt.xlabel('$s$')
        plt.ylabel(r'$E\,\Gamma(s) p(s)$')
        plt.title(r'Low tau, $\beta=1\,\mathrm{meV}^{-1}$')
        plt.xlim([-40,40])
        #plot_cont(f12,6,1,200000,'bdE','etot')
        plt.figure(6)  
        plt.xlabel('$s$')
        plt.ylabel(r'$E\,\Gamma(s) p(s)$')
        plt.title(r'ashoori, $\beta=1\,\mathrm{meV}^{-1}$')
        plt.xlim([-40,20])
        plot_cont(f2,7,1,200000,'bdE','etot')
        
    #plot_s_int(f11,4)
    #plot_s_int(fi,5)
    #plot_s_int(f14,8,opt='log')
    #plot_s_int(f15,9,opt='log')
    """plot_s_int(f15,9)
    plt.ylim([0,10])
    plt.xlim([-5,5])"""
    if 0:
        plot_cont(fi,6,1,200000,'bdE','etot',conn=0)
        plot_cont(fi,7,1,200000,'bdE',conn=0)
        
    
        plt.figure(6)
        plt.xlim([-20,20])
        plt.ylim([0,0.9])
        plt.ylabel(r'$E(s)(1+\mathrm{e}^{-s})p(s)$') 
        plt.title(r'$\mathrm{Without~Metadynamics}$')
        """plt.figure(7)
        plt.xlim([-20,20])
        plt.ylabel(r'$E(s)(1+\mathrm{e}^{-s})p(s)$') 
        plt.title(r'$\mathrm{With~Metadynamics}$')
        plt.ylim([0,0.9])"""
    #plot_cont(f6,7,-1,200000,'bdE')
    #plt.title('New cv')
   #plot_rAB(f0+interac[1]+sym[2]+P[3]+s[0]+'woMetaD/',2,0,'k')

    #plot_energies_vs_t(f6,7,100000)

    
    """plt.figure(0)
    plt.clf()
    Ts = np.linspace(0.1,80,20)
    betas = 1.0/(kB*Ts)
    es = np.zeros(20)
    for i,beta in enumerate(betas):
        es[i] = half_energy(0.1,2,beta,3.0,'fer')
    plt.plot(1.0/(kB*betas),es,'o-')"""