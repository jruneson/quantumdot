import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
from matplotlib import rcParams, interactive
from scipy.special import jv
import copy

#from theoretical_energy import half_energy

plt.rc('text',usetex=True)
plt.rc('font',family='helvetica')
rcParams.update({'font.size': 28})

kB = 1/11.6045
hw = 4.828
#hw=1.21
en_scale = 1#hw/2
offset = 2*hw*0

def plot_cv(f,fig_nr,n=0,opt=None):
    if(n==0):
        data = np.loadtxt(f+'cv.dat')
    else:
        data = np.genfromtxt(f+'cv.dat',max_rows=n)
    #data = np.loadtxt(f+'cv.dat')
    t = data[:,0]
    cv = data[:,1]
    plt.figure(fig_nr)
    plt.clf()
    plt.plot(t,cv,label='s')
    #plt.plot(t,data[:,2])
    plt.xlabel('$t$ (ps)')
    if opt is None:
        plt.ylabel('$s$')
    if(opt=='rw'):
        data = np.loadtxt(f+'rew_factor.dat')
        rw = data[:,1]
        plt.plot(t,10*rw,label='10*rew-factor')
        plt.legend(loc='upper right',fontsize=20)
    if(opt=='bdE'):
        bdE = data[:,2]
        plt.plot(t,bdE,label=r'$\beta\Delta E$')
        plt.legend(loc='upper right',fontsize=20)
        
        
    
def plot_gauss_data(f,fig_nr,opt='Ws',xlim=None,name=None):
    centers = np.loadtxt(f+'cv_centers.dat')
    heights = np.loadtxt(f+'heights.dat')
    plt.figure(fig_nr)
    plt.clf()
    if(opt=='Ws'):
        plt.plot(centers[:,1],heights[:,1],'x')
        plt.xlabel('$s_k$')
        plt.ylabel('$W_k~(k_\mathrm{B}T)$')
    if(opt=='st'):
        plt.plot(centers[:,0]/1000,centers[:,1])
        plt.xlabel('$t\,(ns)$')
        plt.ylabel('$s_k$')
    if(opt=='Wt'):
        plt.plot(heights[:,0]/1000,heights[:,1],'x')
        plt.xlabel('$t\,(ns)$')
        plt.ylabel('$W_k$')
    if(opt=='Wst'):
        plt.plot(centers[:,0]/1000,centers[:,1],'x',label='$s_k$')
        plt.plot(heights[:,0]/1000,100*heights[:,1],'x',label='$100W_k$')
        plt.xlabel('$t\,(ns)$')
        plt.legend(loc='upper right',fontsize=22)
    if xlim is not None:
        plt.xlim(xlim)
    if name is not None:
        plt.title(name)


def plot_cont_sint(f,fig_nr,block_size,opt='Gamma'):
    cv_f = open(f+'cv.dat')
    rw_f = open(f+'rew_factor.dat')
    plt.figure(fig_nr)
    plt.clf()
    plt.xlabel(r'$s$')
    if(opt=='Gamma'):
        plt.ylabel(r'$\Gamma(s)~p(s)$')
    elif(opt=='FES'):
        plt.ylabel(r'$F(s)=-\frac{1}{\beta}\log\,p(s)$')
    interactive(True)
    smin = -60
    smax = 150
    nbins = 500
    beta = 1.0
    cvs = np.linspace(smin,smax,nbins)
    ds = cvs[1]-cvs[0]
    his_rew = np.zeros(nbins)
    #line1, = ax.plot(cvs,his_rew)
    count = 0
    for cv_line, rw_line in zip(cv_f,rw_f):
        if(opt=='Gamma' or opt=='FESbdE'):
            s = float(cv_line.split()[2])
        elif(opt=='FES'):
            s = float(cv_line.split()[1])
        rw = float(rw_line.split()[1])
        b = np.round((s-smin)/(smax-smin)*(nbins-1))
        if b>=nbins:
            b=nbins-1
        if b<0:
            b=0
        his_rew[b] += rw
        count += 1
        if(count==block_size):
            his_rew_norm = his_rew / (sum(his_rew)*ds)
            plt.clf()
            plt.xlabel(r'$s$')
            #plt.plot(cv,his_rew_norm)
            t = cv_line.split()[0]
            if(opt=='Gamma'): 
                plt.plot(cvs,his_rew_norm*(1-np.exp(-cvs)),label=str(t)+' ps')
                plt.ylabel(r'$\Gamma(s)~p(s)$')
            elif(opt=='FES' or opt=='FESbdE'):
                plt.plot(cvs,-np.log(his_rew_norm)/beta,label=str(t)+' ps')                
                plt.ylabel(r'$F(s)=-\frac{1}{\beta}\log\,p(s)$')
            plt.legend(loc='upper right')
            plt.pause(0.01)            
            count = 0
    plt.clf()
    plt.xlabel(r'$s$')    
    #plt.plot(cvs,his_rew_norm)
    if(opt=='Gamma'):
        plt.plot(cvs,his_rew_norm*(1-np.exp(-cvs)))
        plt.ylabel(r'$\Gamma(s)~p(s)$')
        print(sum(his_rew_norm*(1-np.exp(-cvs)))*ds)
    elif(opt=='FES' or opt=='FESbdE'):
        plt.plot(cvs,-np.log(his_rew_norm)/beta,label=r'$-\frac{1}{\beta}\log\,p(s)$')
        plt.ylabel(r'$F(s)~\mathrm{(meV)}$') 
        
def plot_cont(f, fig_nr,sign,block_size,cv='bdE',obs='etot'):
    g_f = open(f+'exc_factor.dat')
    cv_f = open(f+'cv.dat')
    rw_f = open(f+'rew_factor.dat')
    if(obs=='etot'):
        obs_f = open(f+'Total_energy.dat')
    plt.figure(fig_nr)
    plt.clf()
    plt.xlabel(r'$s$')
    plt.ylabel(r'$\Gamma(s)~p(s)$')
    interactive(True)
    if(sign==1):
        if(cv=='G'):
            Gmin = 0
            Gmax = 1000
            file=cv_f
        if(cv=='bdE'):
            Gmin = -50
            Gmax = 40
            file=cv_f
    if(sign==-1):
        if(cv=='G'):
            Gmin = -20
            Gmax = 2
            file = cv_f
        if(cv=='bdE'):
            Gmin = -50
            Gmax = 40
            file= cv_f
    nbins = 1000
    Gs = np.linspace(Gmin,Gmax,nbins)
    dG = Gs[1]-Gs[0]
    his_rew = np.zeros(nbins)
    #line1, = ax.plot(cvs,his_rew)
    count = 0
    for line, rw_line,obs_l in zip(file,rw_f,obs_f):
        if(cv=='G'):
            G = float(line.split()[1])
            G = 1+sign*np.exp(-G)
        elif(cv=='bdE'):
            G = float(line.split()[2])
        rw = float(rw_line.split()[1])
        ob = float(obs_l.split()[1])
        rw *= ob
        b = np.round((G-Gmin)/(Gmax-Gmin)*(nbins-1))
        if b>=nbins:
            b=nbins-1
        if b<0:
            b=0
        t = float(line.split()[0])
        his_rew[b] += rw
        count += 1
        if(count==block_size):
            his_rew_norm = his_rew / (sum(his_rew)*dG)
            plt.clf()
            #plt.plot(cv,his_rew_norm)
            #t = line.split()[0]
            if(cv=='G'):
                plt.plot(Gs,np.log(his_rew_norm+1e-9),label=str(t)+' ps')
                plt.xlabel(r'$\Gamma$')
                plt.ylabel(r'$\log\,p(\Gamma)$')            
            elif(cv=='bdE'):
                plt.plot(Gs,his_rew_norm*(1+sign*np.exp(-Gs)),label=str(t)+' ps')
                plt.xlabel(r'$\beta\Delta E$')    
                plt.ylabel(r'$[1-\exp(-\beta\Delta E)]p(\beta\Delta E)$')
                #plt.ylim([-0.1,0.4])
            if(sign==1):
                plt.legend(loc='upper right')
            if(sign==-1):
                plt.legend(loc='upper left')
            plt.pause(0.01)            
            count = 0
    plt.clf()
    #plt.plot(cvs,his_rew_norm)
    if(cv=='G'):
        plt.plot(Gs,np.log(his_rew_norm+1e-9))
        plt.xlabel(r'$\Gamma$')    
        plt.ylabel(r'$\log\,p(\Gamma)$')
    if(cv=='bdE'):
        plt.plot(Gs,his_rew_norm*(1+sign*np.exp(-Gs)))
        plt.xlabel(r'$\beta\Delta E$')
        plt.ylabel(r'$[1-\exp(-\beta\Delta E)]p(\beta\Delta E)$')
        #plt.ylim([-0.1,0.25])
    print(sum(his_rew_norm*dG))
    plt.grid()
    
    
def read_data(f):
    cv_data = np.loadtxt(f+'cv.dat')
    rf_data = np.loadtxt(f+'rew_factor.dat')
    exc_data = np.loadtxt(f+'exc_factor.dat')
    return cv_data[:,0], cv_data[:,2], rf_data[:,1], exc_data[:,1]
    
def plot_cv_hist(fig_nr,cv,rf,exc):
    smin = min(cv)
    smax = max(cv)
    nbins = 1000
    cvs = np.linspace(smin,smax,nbins)
    #his = np.zeros(nbins)
    his_rew = np.zeros(nbins)
    for i,s in enumerate(cv):
        b = int(np.floor((s-smin)/(smax-smin)*(nbins-1)))
        #his[b] += 1.0
        his_rew[b] += rf[i]
    ds = cvs[1]-cvs[0]
    #his /= (sum(his)*ds)
    his_rew /= (sum(his_rew)*ds)
    #plt.plot(cvs,his*(1-np.exp(-cvs)),label='non-reweighted')
    plt.plot(cvs,his_rew*(1-np.exp(-cvs)),label='$\Gamma(s)p(s)$')
    print(ds*sum(his_rew*(1-np.exp(-cvs))))
    xmax = 80
    xmin = -50
    plt.xlim([xmin,xmax])
    plt.ylim([-0.1,0.3])
    plt.xlabel('$s$')
    plt.ylabel('$\Gamma(s) p(s)$')
    plt.legend(loc='upper right',fontsize=18)
    
def plot_energies(f,fig_nr,clear=1,var='beta',linestyle='-',plot_all=False):
    data = np.loadtxt(f+'results.dat',comments='%')
    if(data.ndim==1):
        data = [data,float('NaN')*np.ones(len(data))]
    if(data.ndim>1):
        x = data[:,0]
        epot = data[:,1]/en_scale
        epot_e = data[:,2]/en_scale
        ekin = data[:,3]/en_scale
        ekin_e = data[:,4]/en_scale
        etot = data[:,5]/en_scale-offset
        etot_e = data[:,6]/en_scale
        evir = data[:,7]/en_scale
        evir_e = data[:,8]/en_scale    
        plt.figure(fig_nr)
        if(clear):
            plt.clf()
        if(var=='P'):
            plt.xlabel(r'$P$')
        elif(var=='beta'):
            plt.xlabel(r'$T\,(\mathrm{K})$')
            x = 11.6045/x
        elif(var=='tau'):
            plt.xlabel(r'$1/\tau~(\mathrm{meV})$')
        plt.errorbar(x,etot,etot_e,marker='^',color='k',label='Total energy',linestyle=linestyle)
        if(plot_all==True):   
            plt.errorbar(x,epot,epot_e,marker='x',color='b',label='Potential energy',linestyle=linestyle)
            plt.errorbar(x,ekin,ekin_e,marker='o',color='r',label='Kinetic energy',linestyle=linestyle)
            plt.errorbar(x,evir,evir_e,marker='v',color='g',label='Virial kinetic',linestyle=linestyle)   

        #plt.ylabel(r'$\mathrm{Energy}~(\hbar\omega_0/2)$')
        plt.ylabel(r'$\mathrm{Energy}~(\mathrm{meV})$')
        #plt.ylim([0.7, 2.6])
        #plt.title('Energy difference-CV, $k_\mathrm{B}T=1\,\mathrm{meV}$')
        #plt.title('Distance-corrected CV, $k_\mathrm{B}T=1\,\mathrm{meV}$')
        return x,etot
    """    
    time = data[0]*0.001
    epot = data[1]/en_scale
    epot_e = data[2]/en_scale
    ekin = data[3]/en_scale
    ekin_e = data[4]/en_scale
    evir = data[7]/en_scale
    evir_e = data[8]/en_scale    
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    if(P!=0):
        time=P
        plt.xlabel('$P$')
    plt.errorbar(time,epot,epot_e,marker='x',color=color,label='Potential energy')
    plt.errorbar(time,ekin,ekin_e,marker='o',color=color,label='Kinetic energy')
    plt.errorbar(time,evir,evir_e,marker='v',color=color,label='Virial energy')   
    if(P==0):   
        plt.xlabel('Time (ns)')
    plt.ylabel('Energy $(\hbar\omega_0/2)$')
    return data[0]*0.001"""

def plot_fes(f,fig_nr):
    centers = np.loadtxt(f+'cv_centers.dat')[:,1]
    heights = np.loadtxt(f+'heights.dat')[:,1]
    plt.figure(fig_nr)
    sigma = 4.0
    gamma = 5.0
    cv = np.linspace(centers.min()-sigma,centers.max()+sigma,200)
    V = np.zeros(len(cv))
    for i,s in enumerate(cv):
        for j in range(len(heights)):
            V[i] += gauss(s,centers[j],sigma,heights[j])
    plt.plot(cv,-gamma/(gamma-1.0)* (V-max(V)),label=r'$-\frac{\gamma}{\gamma-1}V(s)$')
    plt.legend(loc='upper right',fontsize=20)
    
def gauss(x,mu,sigma,h):
    return h*np.exp(-(x-mu)**2/(2.0*sigma**2))

def plot_coord(f,fig_nr):
    data = np.loadtxt(f+'X_coord_n1p1.dat')
    t = data[:,0]
    x = data[:,1]
    plt.figure(fig_nr)
    plt.clf()
    plt.plot(t,x)
    plt.xlabel('$t (ps)$')
    plt.ylabel('$x (a_0)$')
    
def normalize(p,perr,r,d,norm_shell=True):
    dr = r[1]-r[0]
    norm = p.sum()*dr
    p /= norm 
    perr /= norm
    if(norm_shell):
        norm_shell = np.pi*((r+dr)**d-r**d)
        p /= norm_shell
        perr /= norm_shell
    
def plot_rAB_dist(settings,ax=None):
    f = settings[0]
    d = settings[1]
    name = settings[2]
    marker = settings[4]
    color= settings[5]
    #f = f_name[0]
    #name = f_name[1]
    data = np.loadtxt(f+'Prob_distribution.dat')
    r = data[:,0]
    p = data[:,1]
    p_err = data[:,2]
    """p1 = data[:,3]
    p1_err = data[:,4]
    p2 = data[:,5]
    p2_err = data[:,6]"""
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
    #plt.title('Two fermions, Lennard--Jones-potential')
    
def load_lines(filename,num_lines):
    with open(filename) as file:
        head = [next(file) for x in range(num_lines)]
    for line in head:
        myarray = np.fromstring(line)
    print(head)
    return head
    
def plot_energies_vs_t(f,fig_nr,n=50000):
    data = np.genfromtxt(f+'Pot_energy.dat',max_rows=n)
    epot = data[:,1]
    data = np.genfromtxt(f+'Kinetic_energy.dat',max_rows=n)
    ekin = data[:,1]    
    data = np.genfromtxt(f+'Kin_en_virial.dat',max_rows=n)
    evir = data[:,1]
    #data = load_lines(f+'Total_energy_cl.dat',n)
    data = np.genfromtxt(f+'Total_energy.dat',max_rows=n)
    etot = data[:,1]
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
    plt.figure(fig_nr)
    plt.clf()
    plt.plot(t,etot,marker='.',label='etot')
    plt.plot(t,epot,marker='.',label='epot')
    #plt.plot(t,ekin,marker='.',label='ekin')
    plt.plot(t,evir,marker='.',label='evir')
    #plt.plot(t,etotc,marker='.',label='etotc')
    #plt.plot(t,esprc,marker='.')
    #plt.xlim([0,10])
    #plt.plot(t,ekin)
    #plt.plot(t,evir)

def get_axis_limits(ax, scale=0.85):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale

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
    ax1.xaxis.set_major_locator(tmp)
    
def plot_rAB(f,fig_nr,clear=True,color='blue',marker='x',name=None,linestyle='-',show_errors=1):
    if name is None:
        name=f[-25:]
    data = np.loadtxt(f+'Prob_distribution.dat')
    r = data[:,0]
    p = data[:,1]
    p_err = data[:,2]
    normalize(p,p_err,r,2)
    p *= 1e8
    p_err *= 1e8
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    plt.plot(r,p,color=color,linestyle=linestyle)
    if(show_errors==1):
        plt.errorbar(r[::],p[::],p_err[::],linestyle='None',label=name,marker=marker,color=color)
    if(show_errors==3):
        plt.errorbar(r[::3],p[::3],p_err[::3],linestyle='None',label=name,marker=marker,color=color)
    plt.xlabel('$r~(a_0)$')
    plt.ylabel('$p(r_\mathrm{AB})/10^{-8}$')
    return r,p
    
def plot_rAB_th(fig_nr,r,sym):
    P = 10
    hwb=3.0
    hwt = hwb/P
    f = 1+0.5*hwt**2+0.5*hwt*np.sqrt(hwt**2+4)
    c = 1/((f+1)/np.sqrt(f) * (f**P-1)/(f**P+1))
    print(c)
    a = np.sqrt(1.0/(3.675e-5*3.0))
    x = np.linspace(0,10*a,1000)
    dx = x[1]-x[0]
    p2part = copy.deepcopy(r)
    for i,R in enumerate(r):
        xi = -2j*x*R/a**2
        exponential = np.exp(-(2*x**2+R**2)/a**2)
        if(sym=='dis'):
        #p2part = np.exp(-0.5*(r/a)**2)
            bessels = jv(0,xi)+0.5*((2*x**2+R**2)*jv(0,xi)+1j*2*x*R*jv(1,xi))/a**2 * np.exp(-hwb)
            p2part[i]=np.absolute(sum(x*R*exponential*bessels)*dx)
            #p2part[i]=R**2*np.exp(-0.5*R**2/a**2)
            color = 'r'
        if(sym=='bos'):
            p2part[i]=np.absolute(sum(x*R*exponential*jv(0,xi))*dx)
            #p2part[i]=np.absolute(sum(x*R*exponential*(2*np.pi*jv(0,2*xi)+jv(0,xi)**2))*dx)
            color = 'b'
        if(sym=='fer'):
            p2part[i]=np.absolute(sum(x*R**3*exponential*jv(0,xi))*dx)
            color = 'g'
        #p2part = r**2*np.exp(-(r/a)**2)
    normalize(p2part,0,r,2,1)
    scale=1e6
    plt.figure(fig_nr)
    plt.plot(r,p2part*scale,color=color,marker='.',label='Theory Fermion')
    

if __name__=="__main__":
    f0 = '../'
    f1 = '../run1/'
    f2 = '../run2/'
    f3 = '../run3/'
    f4 = '../run4/'
    f5 = '../run5/'
    f6 = '../run6/'
    f7 = '../run7/'
    f8 = '../run8/'
    f9 = '../run9/'
    fs=[f1,f2,f3,f4,f5]
    interac = ['noInt/','LJ/']
    sym = ['Dist/','Boson/','Fermion/']
    P=['P10/','P20/','P40/','P50/','P60/']
    s=['noLJ/','s20/','s50/','s100/']
    dt=['','dt5e-3/','dt1e-3/','dt1e-4/']
    t=['t100000/']
    ff = f0+interac[1]+sym[2]
    f=f0+interac[0]+sym[2]
    flab = '../workstation_lab/metaD_test5/'
    flabb= '../workstation_lab/metaD_test4_20ns/'
    flab1= flab+'run1/'
    flab2= flab+'run2/'
    flab3= flab+'run3/'
    flab4= flab+'run4/'
    flab5= flab+'run5/'
    flab6= flab+'run6/'
    flab7= flab+'run7/'
    flab8= flab+'run8/'
    flab9= flab+'run9/'
    flab10=flab+'run10/'
    
    fi = flab+'run1/'
    if(f!='.'):
        plot_gauss_data(f6,7,'Wt',name='fer')
        plot_gauss_data(f7,8,'Wt',name='bos')
        #plot_rABs(fs)
        plot_energies_vs_t(f6,1,n=50000)
        #plt.legend(loc='upper right',fontsize=18)
        #plot_rAB(f,2,1,'b','^',name='Dist')
        [r,p1]=plot_rAB(f7,2,1,'r','o',name='Singlet')
        [r,p2]=plot_rAB(f6,2,0,'g','v',name='Triplet')
        plt.plot(r,(p1-p2),color='k',marker='.',label='S$-$T')
        plt.ylim([-1,4])
        #plt.grid(True)
        #r=plot_rAB(f4,2,0,'c','v',name='dt=0.1fs')
        #plot_rAB_th(1,r,'fer')
        plt.legend(loc='upper right',fontsize=18)
        plt.title(r'$\mathrm{With MetaD},\quad \beta=1\,\mathrm{meV}^{-1}$')
        #plt.title(r'$\mathrm{No~interaction}$')
        #plot_energies(f1,3,1,'b')
        #plot_energies(f2,3,1,'g')
        #plot_energies(f3,3,1,'r')
        r,ef = plot_energies(f1,3,1,'beta','-',True)
        plt.legend(loc='upper left',fontsize=18)
        r,eb = plot_energies(f2,3,0,'beta','--',True)
        #plt.plot(r,ef-eb,'y')
        """plot_energies(f4,3,0,'r',P=20)
        plot_energies(f5,3,0,'k',P=20)
        plot_energies(f6,3,0,'r',P=50)
        plot_energies(f7,3,0,'g',P=50)
        plot_energies(f8,3,0,'b',P=5)
        plot_energies(f9,3,0,'g',P=30)"""
        #plt.title(r'$\mathrm{Boson}$')
        plt.title(r'$\mathrm{GaAs,}\quad R_W=1.34$')
        plt.xlim([0,40])
        plt.ylim([0,32])
        
        #plot_cv(f3,4,200000)        
        
        if 0:
            plot_cont_sint(f6,3,200000,'FES')
            plt.title('fer')
            plot_cont_sint(f7,4,200000,'FES')
            plt.title('bos')
        if 1:
            #plot_cont(f6,5,-1,200000,'bdE','etot')
            plt.figure(5)
            plt.ylabel(r'$E\, \Gamma(s) p(s)$')
            plt.title('Fermion')
            #plot_cont(f7,6,1,200000,'bdE','etot')
            plt.figure(6)            
            plt.ylabel(r'$E\, \Gamma(s) p(s)$')
            plt.title('Boson')
        if 0:
            plot_cont(f5[0],5,-1,200000,'G')
        #plot_cont(f6,7,-1,200000,'bdE')
        #plt.title('New cv')
        #plot_rAB(f0+interac[1]+sym[2]+P[3]+s[0]+'woMetaD/',2,0,'k')

    plt.show()
    
    """plt.figure(0)
    plt.clf()
    Ts = np.linspace(0.1,80,20)
    betas = 1.0/(kB*Ts)
    es = np.zeros(20)
    for i,beta in enumerate(betas):
        es[i] = half_energy(0.1,2,beta,3.0,'fer')
    plt.plot(1.0/(kB*betas),es,'o-')"""

    
    for i in plt.get_fignums():
        plt.figure(i)
        plt.tight_layout()
