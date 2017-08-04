import numpy as np
from scipy import ndimage,misc,interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as tic
from matplotlib import rcParams, interactive, cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import jv
from scipy.fftpack import fft
from matplotlib.ticker import FormatStrFormatter

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import copy



#from theoretical_energy import half_energy

plt.rc('text',usetex=True)
plt.rc('font',family='helvetica')
rcParams.update({'font.size': 28})

kB = 1/11.6045
hw = 3.0 #4.828
#hw=1.21
en_scale = hw
offset = 2*hw*0

  
def plot_cont(f, fig_nr,sign,block_size,cv='bdE',obs=None,conn=False):
    if(conn):
        conn_sign=1
    else:
        conn_sign=-1
    g_f = open(f+'exc_factor.dat')
    cv_f = open(f+'cv.dat')
    rw_f = open(f+'rew_factor.dat')
    #obs_f = g_f#open(f+'Total_energy.dat')
    plt.figure(fig_nr)
    plt.clf()
    plt.xlabel(r'$s$')
    plt.ylabel(r'$\Gamma(s)~p(s)$')
    interactive(True)
    if(sign==1):
        if(cv=='G'):
            Gmin = 0
            Gmax = 1000
            file=g_f
        if(cv=='bdE'):
            Gmin = -51
            Gmax = 51
            file=cv_f
    if(sign==-1):
        if(cv=='G'):
            Gmin = -500
            Gmax = 1.1
            file = g_f
        if(cv=='bdE'):
            Gmin = -51
            Gmax = 51
            file= cv_f
    nbins = 200
    Gs = np.linspace(Gmin,Gmax,nbins)
    dG = Gs[1]-Gs[0]
    his_rew = np.zeros(nbins)
    #line1, = ax.plot(cvs,his_rew)
    count = 0
    for line, rw_line in zip(file,rw_f):
        if(cv=='G'):
            G = float(line.split()[1])
            #G = 1+sign*np.exp(-G)
        elif(cv=='bdE'):
            G = float(line.split()[2])
        rw = float(rw_line.split()[1])
        if(obs=='etot'):
            ob = float(obs_l.split()[2])
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
                plt.xlabel(r'$W$')
                plt.ylabel(r'$\log\,p(W)$')            
            elif(cv=='bdE'):
                plt.plot(Gs,his_rew_norm*(1+sign*np.exp(conn_sign*Gs)),label=str(t)+' ps')
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
        plt.plot(Gs,np.log(his_rew_norm))
        plt.xlabel(r'$W$')    
        plt.ylabel(r'$\log\,p(W)$')
        plt.ylim([-20,1])
    if(cv=='bdE'):
        plt.plot(Gs,his_rew_norm*(1+sign*np.exp(conn_sign*Gs)))
        plt.xlabel(r'$s$')
        plt.ylabel(r'$(1-\mathrm{e}^{-s})p(s)$')        
        #plt.xlabel(r'$\beta\Delta E$')
        #plt.ylabel(r'$[1-\exp(-\beta\Delta E)]p(\beta\Delta E)$')
        #plt.ylim([-0.1,0.25])
    #print(sum(his_rew_norm*(1+sign*np.exp(conn_sign*Gs))*dG))
    plt.grid()
    return Gs,np.log(his_rew_norm)


def plot_cv(f,fig_nr,n=0,opt=None):
    if(n==0):
        data = np.loadtxt(f+'cv.dat')
    else:
        data = np.genfromtxt(f+'cv.dat',max_rows=n)
    #data = np.loadtxt(f+'cv.dat')
    t = data[:,0]/1000
    cv = data[:,1]
    plt.figure(fig_nr)
    plt.clf()
    plt.plot(t,cv,label='s')
    #plt.plot(t,data[:,2])
    plt.xlabel('$t$ (ns)')
    if opt is None:
        plt.ylabel('$s$')
    if(opt=='rew'):
        data = np.genfromtxt(f+'rew_factor.dat',max_rows=n)
        rw = data[:,1]
        plt.plot(t,5*rw,label='rew-factor')
        plt.legend(loc='upper right',fontsize=20)
    if(opt=='bdE'):
        bdE = data[:,2]
        plt.plot(t,1*bdE,label=r'$\beta\Delta E$')
        plt.legend(loc='upper right',fontsize=20)
    """if(opt=='exp'):
        plt.clf()
        plt.plot(t,np.exp(cv))
        plt.plot(t,data[:,3]*100)"""
    if(opt=='perm'):
        plt.plot(t,data[:,3]*10)
    if(opt=='onlyperm'):
        plt.clf()
        plt.plot(t,data[:,3]*10)
        
        
    
def plot_gauss_data(f,fig_nr,clear=True,opt='Ws',xlim=None,name=None):
    centers = np.loadtxt(f+'cv_centers.dat')
    t = centers[:,0]/1000
    c = centers[:,1]
    if np.size(centers,1)>2:
        W = centers[:,2]
    else:
        heights = np.loadtxt(f+'heights.dat')
        W = heights[:,1]
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    if(opt=='Ws'):
        plt.plot(c,W,'xr')
        plt.xlabel('$s_k$')
        #plt.ylabel('$W_k~(k_\mathrm{B}T)$')
        plt.ylabel(r'$H_k~(\mathrm{meV})$')
    if(opt=='st'):
        plt.plot(t,c,'x')
        plt.xlabel('$t\,(ns)$')
        plt.ylabel('$s_k$')
    if(opt=='Wt'):
        plt.plot(t,W,'x')
        plt.xlabel(r'$t\,(\mathrm{ns})$')
        plt.ylabel(r'$H_k~(\mathrm{meV})$')
    if(opt=='Wst'):
        plt.plot(t,c,'x',label='$s_k$')
        plt.plot(t,100*W,'x',label='$100W_k$')
        plt.xlabel(r'$t\,(\mathrm{ns})$')
        plt.legend(loc='upper right',fontsize=22)
    if xlim is not None:
        plt.xlim(xlim)
    if name is not None:
        plt.title(name)


def plot_cont_sint(f,fig_nr,block_size,beta,opt='Gamma'):
    cv_f = open(f+'cv.dat')
    rw_f = open(f+'rew_factor.dat')
    plt.figure(fig_nr)
    #plt.clf()
    plt.xlabel(r'$s$')
    if(opt=='Gamma'):
        plt.ylabel(r'$\Gamma(s)~p(s)$')
    elif(opt=='FES'):
        plt.ylabel(r'$F(s)=-\frac{1}{\beta}\log\,p(s)$')
    interactive(True)
    smin = -60
    smax = 150
    nbins = 500
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
            """plt.clf()
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
            plt.pause(0.01)"""            
            count = 0
    #plt.clf()
    plt.xlabel(r'$s$')    
    #plt.plot(cvs,his_rew_norm)
    print(count)
    if(opt=='Gamma'):
        plt.plot(cvs,his_rew_norm*(1-np.exp(-cvs)))
        plt.ylabel(r'$\Gamma(s)~p(s)$')
        print(sum(his_rew_norm*(1-np.exp(-cvs)))*ds)
    elif(opt=='FES' or opt=='FESbdE'):
        plt.plot(cvs,-np.log(his_rew_norm)/beta,label=r'$-\frac{1}{\beta}\log\,p(s)$')
        plt.ylabel(r'$F(s)~\mathrm{(meV)}$') 
        return cvs,-np.log(his_rew_norm)/(beta)


def bin_hist(x,y):
    x = x[:-1:2]
    y = 0.5*(y[:-1:2]+y[1::2])
    return x,y
        
def plot_s_int(f, fig_nr, clear=1, opt='',color='b',label='',with_wall=False,linestyle='-',linewidth=1,ax=None,marker='',markersize=1,withMetaD=True):
    data = np.loadtxt(f+'CV_distributions.dat')
    s = data[:,0]
    shist = data[:,1]
    Whist = data[:,2]
    Ehist = data[:,3]
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    #for i in range(2):
    #    s,hist=bin_hist(s,shist)
    if(opt=='log'):
        linewidth=1
        if not with_wall:
            linewidth=1.5
        plt.plot(s,-np.log(shist))
        #plt.plot(s,-np.log(shist),label=labels[0],linestyle=linestyle,linewidth=linewidth)
        if data.shape[1]>4:        
            bias = data[:,4]
            if with_wall:
                wall = np.zeros(len(s))
                for i,S in enumerate(s):
                    if S>30:
                        wall[i] = 0.2/2*(S-30)**2
            V = -bias+np.max(bias)
            #plt.plot(s,V,label='-V(s)')
            if with_wall:
                plt.plot(s,-np.log(shist)+wall,label=labels[1],linestyle='--')
            #plt.plot(s,-np.log(shist)+bias-np.max(bias),label='F(s)+V(s)')
            if data.shape[1]>5 and 0:
                bias_der = data[:,5]
                plt.plot(s,bias_der,label='grad V')
            plt.legend(loc='upper center',fontsize=20)
    elif(opt=='shist'):
        plt.plot(s,shist,color=color,label=label,linewidth=linewidth)
    elif(opt=='Whist'):
        if ax is None:
            plt.plot(s,Whist,color=color,label=label)
        else:
            if(withMetaD ):
                maxi=1025
                mini=975
                maxi2=1010
                mini2=990
                ax.plot(s[mini2:maxi2],Whist[mini2:maxi2],color=color,label=label,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize)
                ax.plot(s[mini:mini2:4],Whist[mini:mini2:4],color=color,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize)
                ax.plot(s[maxi2:maxi:4],Whist[maxi2:maxi:4],color=color,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize)
                ax.plot(s[:mini:10],Whist[:mini:10],color=color,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize)
                ax.plot(s[maxi::10],Whist[maxi::10],color=color,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize)
            else:
                ax.plot(s,Whist,color=color,label=label,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize)

        return s,Whist
    else:
        plt.plot(s,shist,'b:', label=r'$p(s)$',linewidth=2)
        plt.plot(s,Whist,'g--', label=r'$W p(s)$')
        plt.plot(s,Ehist,'r-',label=r'$E W p(s)$')
    plt.xlim([-50,50])
    
        
    
    
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
    
def plot_energies_vs_t(f,fig_nr,n=100000,P=20,conn=False):
    data = np.genfromtxt(f+'Total_energy.dat',max_rows=n)
    signal = data[:,2]
    t = data[:,0]
    plt.figure(fig_nr)
    plt.clf()
    plt.plot(t,signal)
    """dt = t[1]-t[0]
    N = len(signal)
    xf = np.linspace(0.0,1.0/(2.0*dt),N/2)
    yf = np.fft.fft(signal)*2.0/N
    print(N)
    plt.plot(xf,np.abs(yf[:N//2]))
    plt.xlim([0,15])
    plt.ylim([0,6])
    plt.xlabel(r'$f~(\mathrm{ps}^{-1})$')
    plt.ylabel(r'$\mathrm{FFT~of~} E$')
    plt.ylabel(r'$\mathrm{FFT~of~}\langle WE\rangle /\langle W \rangle$')"""
    
def plot_energies(f,fig_nr,clear=1,var='beta',beta=1,linestyle='',marker='v',label='Total energy',
                  plot_all=False,col=0,bennett=True,ediff=True,color='b',ax=None):
    if(col==0):
        data = np.loadtxt(f+'results.dat',comments='%')
    elif(bennett):
        data = np.loadtxt(f+'energies_bennett.dat',comments='%')
    else:
        data = np.loadtxt(f+'energies_fermion.dat',comments='%')
    if(data.ndim==1):
        data = [data,float('NaN')*np.ones(len(data))]
    if(data.ndim>1):
        x = data[:,0]
        num_obs = round((np.size(data,1)-1)/2)
        if(num_obs==1):
            etot = data[:,1]/en_scale
            etot_e = data[:,2]/en_scale
        if(num_obs==3):
            etot = data[:,3]/en_scale
            etot_e = data[:,4]/en_scale
        if(num_obs>3):
            etot = data[:,5]/en_scale
            etot_e = data[:,6]/en_scale
            
        if(col>0):
            etot = data[:,col]/en_scale
            etot_e = data[:,col+1]/en_scale
            if(ediff):
                etot += data[:,col+2]/en_scale
                etot_e += data[:,col+3]/en_scale
   
        plt.figure(fig_nr)
        if(clear):
            plt.clf()
        if(var=='P'):
            plt.xlabel(r'$P$')
        elif(var=='beta'):
            #plt.xlabel(r'$T\,(\mathrm{K})$')
            x = 11.6045/x
        elif(var=='tau'):
            x /= beta
            #plt.grid(True)
            plt.xlabel(r'$1/\tau~(\mathrm{meV})$')
        elif(var=='tau2'):
            x = 2.0/x**2
            plt.xlabel(r'$1/\tau^2~(\mathrm{meV}^2)$')
        if ax is None:
            plt.errorbar(x,etot,etot_e,marker=marker,markersize=12,label=label,linestyle=linestyle,color=color)
        else:
            ax.errorbar(x,etot,etot_e,marker=marker,markersize=10,label=label,linestyle=linestyle,color=color)
        if(plot_all==True):   
            epot = data[:,1]/en_scale
            epot_e = data[:,2]/en_scale
            ekin = data[:,3]/en_scale
            ekin_e = data[:,4]/en_scale
            evir = data[:,7]/en_scale
            evir_e = data[:,8]/en_scale 
            plt.errorbar(x,epot,epot_e,marker='x',color='b',linestyle=linestyle)
            plt.errorbar(x,ekin,ekin_e,marker='o',color='r',linestyle=linestyle)
            plt.errorbar(x,evir,evir_e,marker='v',color='g',linestyle=linestyle)   
        #plt.ylabel(r'$E/\hbar\omega_0$')
        #plt.ylabel(r'$\mathrm{Energy}~(\hbar\omega_0)$')
        #plt.ylabel(r'$\mathrm{Energy}~(\mathrm{meV})$')
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
    plt.clf()
    sigma = 2.0
    gamma = 4.0
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
    
    

    
def autocorr(x):
    x -= np.mean(x)
    result = np.correlate(x,x,mode='full')
    return result[result.size//2:]/result[result.size//2]
    


def plot_autocorr(f,fig_nr,clear=True, n=100000,label=''):
    data = np.genfromtxt(f+'Total_energy.dat',max_rows=n)
    t = data[:,0]  - data[0,0] 
    obs = data[:,1]
    corr_f = autocorr(obs)
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    tmax = t[-1]
    plt.plot(t,corr_f,label=label)
    plt.xlim([0,10])
    plt.xlabel(r'$t~(\mathrm{ps})$')
    plt.ylabel(r'$c(t)$')
    plt.legend(loc='upper right',fontsize=24)
    plt.title(r'$\mathrm{Autocorrelation~function}$')

def get_axis_limits(ax, scale=0.85):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale


    
def normalize(p,perr,r,d=2,norm_shell=False):
    dr = r[1]-r[0]
    norm = p.sum()*dr
    p2 = p/norm
    perr2 = perr/norm
    if(norm_shell):
        norm_shell = np.pi*((r+0.01)**d-r**d)
        p2 /= norm_shell
        perr2 /= norm_shell
    return (p2,perr2)
    
    
def plot_rAB(f,fig_nr,clear=True,d=2,color='blue',marker='x',name=None,linestyle='-',show_errors=1,scale=1.0):
    if name is None:
        name=f[-25:]
    data = np.loadtxt(f+'Pair_correlation.dat')
    r = data[:,0]
    p = data[:,1]
    p_err = data[:,2]
    (p2,perr2) = normalize(p,p_err,r,d,1)
    #p *= 1e8
    #p_err *= 1e8
    if(clear):
        plt.clf()
    plt.plot(r,p2*scale,color=color,linestyle=linestyle)
    if(show_errors):
        n = show_errors
        plt.errorbar(r[1::n],p2[1::n]*scale,perr2[1::n]*scale,markersize=10,linestyle='None',label=name,marker=marker,color=color)
    #if(show_errors>1):
    #    n=show_errors
    #    plt.errorbar(r[::3],p[::3],p_err[::3],linestyle='None',label=name,marker=marker,color=color)
    #plt.ylim([-0.1,0.7])
    return r,p2,perr2
    
def plot_rAB_th(fig_nr,r,d,sym):
    m_hbar2 = 0.01323
    hw = 3.0
    beta = 1
    hwb = hw*beta
    a = np.sqrt(1.0/(m_hbar2*hw))
    p2part = copy.deepcopy(r)
    if(d==1):
        if(sym=='dis'):
            #print(np.exp(-hwb))
            prob_exc =  np.exp(-hwb)
            p2part=np.exp(-0.5*(r/a)**2)*(1+(1+0.5*(r/a)**2)*prob_exc)
            label=''
            color='r'
        if(sym=='bos'):
            p2part=np.exp(-r**2/(2*a**2))
            label=''
            color='b'
        if(sym=='fer'):
            p2part=r**2*np.exp(-r**2/(2*a**2))
            label=''
            color='g'
        #color='k'
    else:
        x = np.linspace(0,10*a,1000)
        dx = x[1]-x[0]
        for i,R in enumerate(r):
            xi = 2j*x*R/a**2
            exponential = np.exp(-(2*x**2+R**2)/a**2)
            if(sym=='dis'):
            #p2part = np.exp(-0.5*(r/a)**2)
                bessels = jv(0,xi)+0.5*((2*x**2+R**2)*jv(0,xi) - 2*1j*x*R*jv(1,xi))/a**2 * np.exp(-hwb)
                p2part[i]=(sum(x*R*exponential*bessels)*dx)
                #p2part[i]=R**2*np.exp(-0.5*R**2/a**2)
                label=''
                color='r'
            if(sym=='bos'):
                p2part[i]=(sum(x*R*exponential*jv(0,xi))*dx)
                #p2part[i]=np.absolute(sum(x*R*exponential*(2*np.pi*jv(0,2*xi)+jv(0,xi)**2))*dx)
                label=''
                color = 'b'
            if(sym=='fer'):
                p2part[i]=(sum(x*R**3*exponential*jv(0,xi))*dx)
                label=''
                color='g'
        #color='k'
            #p2part = r**2*np.exp(-(r/a)**2)
    dr = r[1]-r[0]
    #print('dr='+str(dr))
    (p2,_)=normalize(p2part,p2part,r,d,1)
    #print(p2.sum()*dr)
    #scale=1e6
    if(fig_nr>-1):
        plt.figure(fig_nr)
    if(sym=='bos' and d==2):
        plt.plot(r[3:],p2[3:],color=color,label=label)
    else:       
        plt.plot(r[1:],p2[1:],color=color,label=label)
    return p2
    

    
def bias(s,h_arr,cvc_arr):
    value = 0
    for h,cvc in zip(h_arr,cvc_arr):
        value += gaussian(s,h,cvc)
    return value

def gaussian(s,h,cvc):
    sigma = 4
    #print(h)
    #print(np.exp(-(s-cvc)**2/(2*sigma**2)))
    return h*np.exp(-(s-cvc)**2/(2*sigma**2))
            
def f_FD(x):
    return 1.0/(1+np.exp(x))            
            
def plot_1d_dist(fs,fig_nr,clear=1,n=2,ax=None):
    colors = ['b','g','r']
    labels = ['Boson','Fermion','Partial']
    markers= ['s','D','o']
    plt.figure(fig_nr)
    if(clear):
        plt.clf()    
    for i in range(3):
        if i==0:
            f = fs[0]
        else:
            f = fs[1]
        data = np.loadtxt(f+'Prob_dist1d.dat')
        r = data[:,0]
        dim = round((np.size(data,1)-1)/2)

        for d in range(dim):
            p_pre = data[:,2*d+1]
            perr_pre = data[:,2*d+2]
            p,perr = normalize(p_pre,perr_pre,r,d=1)
            if i==0:
                p0 = copy.deepcopy(p)
            if i==2:
                p = p - p0*0.5
            if ax is None:
                plt.errorbar(r[::n],p[::n],perr[::n],markersize=10,color=colors[i],label=labels[i],marker=markers[i])
            else:
                ax.errorbar(r[::n],p[::n],perr[::n],markersize=10,color=colors[i],label=labels[i],marker=markers[i])
    plt.xlim([-15,15])
    plt.ylim([-0.0,0.12])
    return r
    
def plot_1d_dist_th(fig_nr,r,d=1,hw=3.0,ax=None):
    m_hbar2 = 0.013234
    c = m_hbar2*hw
    p0_pre = np.exp(-c*r**2)
    p1_pre = r**2*np.exp(-c*r**2)
    (p0,_) = normalize(p0_pre,p0_pre,r,d)
    (p1,_) = normalize(p1_pre,p1_pre,r,d)
    plt.figure(fig_nr)
    lw = 2
    if ax is None:
        plt.plot(r,p0,'k--')
        plt.plot(r,0.5*(p0+p1),'k--')
        plt.plot(r,0.5*p1,'k--')
    else:
        ax.plot(r,p0,'k--',linewidth=lw)
        ax.plot(r,0.5*(p0+p1),'k--',linewidth=lw)
        ax.plot(r,0.5*p1,'k--',linewidth=lw)

    
def smooth(data, n):
    for k in range(n):
        for i in range(1,data.shape[0]-2):
            for j in range(1,data.shape[1]-2):
                data[i,j] = (data[i-1,j]+data[i,j-1]+4*data[i,j]+data[i+1,j]+data[i,j+1])*0.125
    return data

#def darken(x, ):
#    return x*0.8   
    
def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)
        
        
def plot_2d_dist_temps(f,fig_nr,num_smooths=0):   
    fig = plt.figure(num=fig_nr,figsize=(15,5.7))   
    gs = mpl.gridspec.GridSpec(1,3)
    gs.update(wspace=0.02,bottom=0.16,top=0.85)
    #fig.set_size_inches(12,6.2)
    plt.clf()  
    betas=['beta0-34/','beta0-67/','beta1/']
    titles = [r'$T=34\,\mathrm{K}$',r'$T=17\,\mathrm{K}$',r'$T=12\,\mathrm{K}$']
    for i in range(3):
        woM = ''
        bu = ''
        if i==2:
            woM = 'MetaD/'
            bu = 'backup/'
        data_s = np.loadtxt(f+'singlet/'+bu+betas[i]+woM+'disconnected/Prob_dist2d.dat')
        data_t = np.loadtxt(f+'triplet/'+bu+betas[i]+woM+'Prob_dist2d.dat')
        r1 = data_s[1:,0]
        r2 = data_s[0,1:]
        dr1 = r1[1]-r1[0]
        dr2 = r2[1]-r2[0]
        num_smooths=0
        hist_s = smooth(data_s[1:,1:],num_smooths)
        hist_t = smooth(data_t[1:,1:],num_smooths)
        """for k in range(0,hist_t.shape[0]-1):
            for l in range(0,hist_t.shape[1]-1):
                if hist_t[k,l]<0:
                    hist_t[k,l]=0"""
        norm_s = sum(sum(hist_s))*dr1*dr2
        norm_t = sum(sum(hist_t))*dr1*dr2
        hist_s = hist_s/norm_s*10000
        hist_t = hist_t/norm_t*10000
        X,Y = np.meshgrid(r1,r2)
        if(i==2 and 0):
            hist_t = ndimage.interpolation.zoom(hist_t,2)
        Z_s = hist_s.reshape(X.shape)
        Z_t = hist_t.reshape(X.shape)
        if i==1:
            Z_s = Z_s.transpose()
        if i==2:
            Z_s = Z_s.transpose()
            Z_t = Z_t.transpose()
        Z = Z_s #- Z_s*0.5
        if 1:
            #Z = ndimage.interpolation.zoom(Z,2)
            #Z = misc.imresize(Z,5.0,interp='bilinear')
            Z = ndimage.gaussian_filter(Z,4)
            #Z = misc.imresize(Z,0.2,interp='bilinear')
            #Z = ndimage.interpolation.zoom(Z,0.5)
        if i==0:
            ax0 = fig.add_subplot(gs[0])
            plt.ylabel(r'$y~(\mathrm{nm})$')
            ax0.set_title(titles[i],fontsize=28)
            Z *= 1.2
        if i==1:
            ax1 = fig.add_subplot(gs[1],sharey=ax0)
            plt.setp(ax1.get_yticklabels(),visible=False)
            ax1.set_title(titles[i],fontsize=28)
            Z *= 1
            
        if i==2:
            ax2 = fig.add_subplot(gs[2],sharey=ax0)
            plt.setp(ax2.get_yticklabels(),visible=False)
            ax2.set_title(titles[i],fontsize=28)
        im = plt.imshow(Z,interpolation='bilinear',extent=[r1[0],r1[-1],r2[0],r2[-1]],
                        cmap='jet',vmin=0,vmax=9)
        plt.ylim([-40,40])
        plt.xlabel(r'$x~(\mathrm{nm})$')
    if 1:
        ax0.set_xticks([-30,-15,0,15,30])
        ax1.set_xticks([-30,-15,0,15,30])
        ax2.set_xticks([-30,-15,0,15,30])
    cax = fig.add_axes([0.905,0.16,0.03,0.69])
    fig.colorbar(im,cax=cax,orientation='vertical',ticks=None)
    #plt.suptitle(r'$\mathrm{Singlet~at~ellipticity}~\omega_y/\omega_x=1.38$',x=0.53,y=0.995,fontsize=32)

        
def plot_2d_dist(folders,fig_nr,titles,suptitle='', stride=1,use_contour=True,
                 num_smooths=0):
    fig = plt.figure(num=fig_nr,figsize=(15,5.7))   
    gs = mpl.gridspec.GridSpec(1,3)
    gs.update(wspace=0.02,bottom=0.16,top=0.85)
    #fig.set_size_inches(12,6.2)
    plt.clf()  
    for i in range(3):
        if i==0:
            f = folders[0]
        else:
            f = folders[1]
        data = np.loadtxt(f+'Prob_dist2d.dat')
        r1 = data[1:,0]
        r2 = data[0,1:]
        dr1 = r1[1]-r1[0]
        dr2 = r2[1]-r2[0]
        histpre = data[1:,1:]
        for k in range(0,histpre.shape[0]-1):
            for l in range(0,histpre.shape[1]-1):
                if histpre[k,l]<0:
                    histpre[k,l]=0
        num_smooths=0
        hist = smooth(histpre,num_smooths)
        norm = sum(sum(hist))*dr1*dr2
        scale = 10000
        hist = hist/norm*scale
        X,Y = np.meshgrid(r1,r2)
        use_filter = False
        Z = hist.reshape(X.shape)
        if 0:
            if i==0:
                Z = Z.transpose()
            if i>=1:
                Z = ndimage.interpolation.zoom(Z,2.0)
        if 1:
            #Z = ndimage.interpolation.zoom(Z,3.0)
            #Z = misc.imresize(Z,5.0,interp='bilinear')
            Z = ndimage.gaussian_filter(Z,3)
            #Z = misc.imresize(Z,0.2,interp='bilinear')
            #Z = ndimage.interpolation.zoom(Z,0.25)
        if i==0:
            ax0 = fig.add_subplot(gs[0])
            #Z = ndimage.interpolation.zoom(Z,5.0)
            #Z = misc.imresize(Z,5.0,interp='bilinear')
            Z0 = copy.deepcopy(Z)
            Z = Z*0.7
        if i==1:
            ax1 = fig.add_subplot(gs[1],sharey=ax0)
            if(use_filter and 0):
                Z = ndimage.gaussian_filter(Z,3)
                #Z = misc.imresize(Z,0.2,interp='bilinear')
            if use_contour:
                plt.setp(ax1.get_yticklabels(),visible=False)
        if i==2:
            ax2 = fig.add_subplot(gs[2],sharey=ax0)
            if(use_filter):
                Z = ndimage.gaussian_filter(Z,3)
            #Z = misc.imresize(Z,0.2,interp='bilinear')
            Z = Z-0.5*Z0
            plt.setp(ax2.get_yticklabels(),visible=False)
            Z *= 2.0
        if use_contour:
                #xticks[0].label1.set_visible(False)
            #levels=np.arange(2,11,2)
            plt.set_cmap('hot')
            #dark_inferno=cmap_map(darken,plt.get_cmap('inferno_r'))
            im = plt.imshow(Z,interpolation='bilinear',extent=[r1[0],r1[-1],r2[0],r2[-1]],
                           vmin=0,vmax=8,cmap=('hot'))
            #c = plt.contour(X,Y,Z,cmap=cm.viridis)#,levels=levels)
            #plt.clabel(c, inline=1,fontsize=16,fmt='%1.1f')
            plt.xlabel(r'$x~(\mathrm{nm})$')
            if i==0:
                plt.ylabel(r'$y~(\mathrm{nm})$')
                ax0.set_title(titles[0],fontsize=28)
            if i==1:
                ax1.set_title(titles[1],fontsize=28)
            if i==2:
                ax2.set_title(titles[2],fontsize=28)
            plt.ylim([-15,15])
            plt.ylim([-40,40])

        else:
            ax = fig.add_subplot(111,projection='3d')
            surf = ax.plot_surface(X,Y,Z, rstride=stride, cstride=stride, cmap=cm.viridis)
            ax.set_zlabel(r'$p(x,y)$')
            ax.set_xlabel(r'$x~(\mathrm{nm})$')
            ax.set_ylabel(r'$y~(\mathrm{nm})$')
            ax.set_zlim([-0.05,0.5])
    #plt.set_aspect('equal')
    cax = fig.add_axes([0.905,0.16,0.03,0.69])
    fig.colorbar(im,cax=cax,orientation='vertical',ticks=None)
    plt.subplots_adjust(wspace=None)
    plt.suptitle(suptitle,x=0.53,y=0.99,fontsize=32)
    if 1:
        #xticksc = ax2.get_xticks()
        #print(xticksc)
        ax0.set_xticks([-30,-15,0,15,30])
        ax1.set_xticks([-30,-15,0,15,30])
        ax2.set_xticks([-30,-15,0,15,30])
    if use_contour:
        ax0.set_aspect('equal')
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
    #fig.tight_layout()
        
def plot_2dpaircorr(folders,fig_nr,titles,suptitle,num_smooths=0,m=1):    
    fig = plt.figure(num=fig_nr,figsize=(11,5.7))   
    gs = mpl.gridspec.GridSpec(1,2)
    gs.update(wspace=0.3,bottom=0.2,top=0.85)
    plt.clf()
    for i in range(2):
        f = folders[i]
        data = np.loadtxt(f+'Pair_corr2d.dat')
        r1 = data[:,0]
        r2 = data[:,1]
        n = np.sqrt(np.size(r1))
        dr1 = r1[n]-r1[0]
        dr2 = r2[1]-r2[0]
        p = data[:,2]
        #X = r1.reshape((n,n))
        Zpre = p.reshape((n,n))
        Z = smooth(Zpre,num_smooths) 
        norm = sum(sum(Z))*dr1*dr2
        print(norm)
        Z = Z/norm
        Z = Z*10000
        if 1:
            Z = Z.transpose()
        if i==0:
            ax0 = fig.add_subplot(gs[0])
            plt.ylabel(r'$y~(\mathrm{nm})$')
            Z0 = Z
            Zx0,Zy0 = twoDto1Dproj(Z,dr1)       
        if i==1:
            ax1 = fig.add_subplot(gs[1],sharey=ax0)
            plt.setp(ax1.get_yticklabels(),visible=False)
            Zx1,Zy1 = twoDto1Dproj(Z,dr1)
        im = plt.imshow(Z,interpolation='gaussian',extent=[r1[0],r1[-1],r2[0],r2[-1]],
                    cmap='hot')
        plt.ylim([-60,60])
        plt.xlabel(r'$x~(\mathrm{nm})$')
        if i==0:
            ax0.set_title(titles[0],fontsize=28)
            cax = fig.add_axes([0.46,0.2,0.03,0.65])
            fig.colorbar(im,cax,cax,orientation='vertical')
        if i==1:
            ax1.set_title(titles[1],fontsize=28)
            cax = fig.add_axes([0.90,0.2,0.03,0.65])
            fig.colorbar(im,cax=cax,orientation='vertical',ticks=None)
            #cax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #plt.subplots_adjust(wspace=None)
    if 1:
        ax0.set_xticks([-75,-50,-25,0,25,50,75])
        ax1.set_xticks([-75,-50,-25,0,25,50,75])
    plt.suptitle(suptitle,x=0.53,y=0.985,fontsize=32)
    ax0.set_aspect('equal')
    ax1.set_aspect('equal')
    plt.figure(fig_nr+1)
    plt.clf()
    r1half = r1[n*n/2:n*n:n]
    plt.plot(r1half[::m],Zx0[::m],'sb-',label='Singlet x',markersize=10)
    plt.plot(r1half[::m],Zx1[::m],'Dg-',label='Triplet x',markersize=10)
    plt.plot(r1half[::m],Zy0[::m],'^r-',label='Singlet y',markersize=10)
    plt.plot(r1half[::m],Zy1[::m],'o-',color='darkorange',label='Triplet y',markersize=10)
    plt.xlim([0,95])
    plt.ylim([0,2.8])
    plt.xlabel(r'$r_i~(\mathrm{nm})$')
    plt.ylabel(r'$g(r_i)~(10^{-2}\,\mathrm{nm}^{-1})$')
    plt.legend(loc='upper right',fontsize=22)
    #plt.title(suptitle,fontsize=30)
    #plt.title(r'$\mathrm{Pair~correlation},~\omega_y/\omega_x=3.0$')
    fig = plt.gcf()
    fig.tight_layout()
    
    
def twoDto1Dproj(Z,dr):
    Zx = np.sum(Z,axis=0)
    Zy = np.sum(Z,axis=1)
    n = len(Zx)
    m = len(Zy)
    Zxmirr = (Zx[0:n/2:1][::-1]+Zx[n/2:n])/2
    Zymirr = (Zy[0:m/2:1][::-1]+Zy[m/2:n])/2
    return Zxmirr*dr/100,Zymirr*dr/100
    
def free_energy_diff_mass(f1,f2,f3,f4,fig_nr,beta=1.0):
    c_hist_oo_m1 = np.loadtxt(f1+'fsum_mass.dat')[:,1:]
    c_hist_oo_m2 = np.loadtxt(f2+'fsum_mass.dat')[:,1:]
    c_hist_O_m1 = np.loadtxt(f3+'fsum_mass.dat')[:,1:]
    c_hist_O_m2 = np.loadtxt(f4+'fsum_mass.dat')[:,1:]
    dE_oo_m1 = np.loadtxt(f1+'Mass_diff_distributions.dat')
    dE_oo_m2 = np.loadtxt(f2+'Mass_diff_distributions.dat')
    dE_O_m1 = np.loadtxt(f3+'Mass_diff_distributions.dat')
    dE_O_m2 = np.loadtxt(f4+'Mass_diff_distributions.dat')

    plt.figure(fig_nr)
    plt.clf()
    s_oo1 = dE_oo_m1[:,0]
    s_oo2 = dE_oo_m2[:,0]
    s_O1 = dE_O_m1[:,0]
    s_O2 = dE_O_m2[:,0]
    ds = s_oo1[1]-s_oo2[0]
    hist_oo_m1 = dE_oo_m1[:,1]
    hist_oo_m1 /= sum(hist_oo_m1)*ds
    hist_oo_m2 = dE_oo_m2[:,1]
    hist_oo_m2 /= sum(hist_oo_m2)*ds
    hist_O_m1 = dE_O_m1[:,1]
    hist_O_m1 /= sum(hist_O_m1)*ds
    hist_O_m2 = dE_O_m2[:,1]
    hist_O_m2 /= sum(hist_O_m2)*ds
    plt.plot(s_oo1, (hist_oo_m1), 'b', label='oo m1')
    plt.plot(s_oo2, (hist_oo_m2), 'r', label='oo m2')
    plt.plot(s_O1, (hist_O_m1), 'g', label='O m1')
    plt.plot(s_O2, (hist_O_m2), 'c', label='O m2')
    plt.legend(loc='lower right',fontsize=24)
    
    cs = c_hist_oo_m1[0,:]    
    plt.figure(fig_nr+1)
    plt.clf()    
    numblocks = len(c_hist_oo_m1[1:,0])
    Czeros = np.zeros(numblocks)   
    Czeros2 = np.zeros(numblocks)   
    for i in range(1,numblocks+1):
        b1u = c_hist_oo_m1[i,:]
        b1c = c_hist_oo_m2[i,:]
        b2u = c_hist_O_m1[i,:]
        b2c = c_hist_O_m2[i,:]
        diffl = np.log(b1u)-np.log(b1c)
        diffl2 = np.log(b2u)-np.log(b2c)
        poly = np.polyfit(cs,diffl,1)
        poly2 = np.polyfit(cs,diffl2,1)
        C = -poly[1]/poly[0]
        C2 = -poly2[1]/poly2[0]
        if 1:        
            plt.plot(cs,diffl,'c')
            plt.plot(cs,diffl2,'m')
            plt.plot(cs,np.polyval(poly,cs),'x--')
            plt.plot(cs,np.polyval(poly2,cs),'x--')
            plt.plot([C],[0],'bo')
            plt.plot([C2],[0],'ro')
        Czeros[i-1]=C        
        Czeros2[i-1]=C2
    dF = Czeros/beta
    dF2 = Czeros2/beta
    dF_t = dF2-dF
    print("dF_oo="+str(np.mean(dF))+"+/-"+str(np.std(dF)))
    print("dF_O="+str(np.mean(dF2))+"+/-"+str(np.std(dF2)))
    print("dF_masses="+str(np.mean(dF_t))+"+/-"+str(np.std(dF_t)))
    
    """num=0
    den=0
    num2=0
    den2=0
    #avg=0
    C = -2.7
    C2 = -3.5
    for i,S in enumerate(s):
        num += f_FD(S+C)*hist_oo_m1[i]
        den += f_FD(-S-C)*hist_oo_m2[i]
        num2 += f_FD(S+C)*hist_O_m1[i]
        den2 += f_FD(-S-C)*hist_O_m2[i]
        #avg += np.exp(-S)*hist_O_m1[i]
    dF = -np.log(num/den)/beta + C
    dF2 = -np.log(num2/den2)/beta + C2
    print(dF)
    print(dF2)"""
            
            
def free_energy_diff(f1,f2,fig_nr,beta=1.0,P=10,opt='FB',deg=1,C=-0.04):
    if opt=='mass':
        c_hist_unconn = np.loadtxt(f1+'fsum_mass.dat')[:,1:]
        c_hist_conn = np.loadtxt(f2+'fsum_mass.dat')[:,1:]        
        #dE_unconn = np.loadtxt(f1+'Mass_diff_distributions.dat')[:,0:2]
        #dE_conn = np.loadtxt(f2+'Mass_diff_distributions.dat')[:,0:2]
    else:
        c_hist_unconn = np.loadtxt(f1+'fsum_N2_P'+str(P)+'.dat')[:,1:]
        c_hist_conn = np.loadtxt(f2+'fsum_N1_P'+str(P)+'.dat')[:,1:]
    dE_unconn = np.loadtxt(f1+'CV_distributions.dat')[:,0:2]
    dE_conn = np.loadtxt(f2+'CV_distributions.dat')[:,0:2]
    #print(c_hist_conn)
    plt.figure(fig_nr)
    plt.clf()
    s = dE_unconn[:,0]
    ds = s[1]-s[0]
    hist_oo = dE_unconn[:,1]
    hist_oo /= sum(hist_oo)*ds
    s2 = dE_conn[:,0]
    hist_O = dE_conn[:,1]
    hist_O /= sum(hist_O)*ds
    plt.plot(s, -np.log(hist_oo), 'b', label='Disc.')
    plt.plot(s2, -np.log(hist_O),'k--',label='Conn.' )
    #C = -0.041#-0.08
    plt.plot(s, hist_oo*f_FD(s+C)*30)
    plt.plot(s2, hist_O*f_FD(-s2-C)*30)
    plt.xlim([-20,20])
    plt.ylim([0,16])
    plt.xlabel(r'$s$')
    plt.ylabel(r'$F(s) = - \log\,p(s)$')
    plt.tight_layout()
    plt.legend(loc='lower right',fontsize=22)
    num=0
    den=0
    sq_oo=0
    sq_O=0
    for i,S in enumerate(s):
        num += f_FD(S+C)*hist_oo[i]
    for i,S in enumerate(s2):
        den += f_FD(-S-C)*hist_O[i]
        #sq_oo += f_FD(S+C)**2*hist_oo[i]
        #sq_O += f_FD(-S-C)**2*hist_O[i]
    #print(num/den)
    dF =  np.log(num/den) + C
    ds = s[1]-s[0]
    avg_oo = num*ds
    avg_O = den*ds
    #print(avg_oo/avg_O)
    sq_oo *= ds
    sq_O *= ds
    #n=200000
    #err = ((sq_O - avg_O**2)/avg_O**2 + (sq_oo-avg_oo**2)/avg_O**2)/n

    #dF_p = -np.log(ds*sum(hist_oo*np.exp(-s)))
    #print("F_O-F_oo, only oo:\t"+str(dF_p))
    dF_FB = -np.log((1-np.exp(dF))/(1+np.exp(dF)))
    #dF_FBp = -np.log((1-np.exp(-dF_p))/(1+np.exp(-dF)))
    #print("F_F-F_B, only oo:\t"+str(dF_FBp))
    #print("Error in DeltaF:\t"+str(np.sqrt(err)))

    print("F_O-F_oo with f_FD:\t"+str(-dF/beta))
    #print("F_F-F_B with f_FD:\t"+str(dF_FB/beta))
    
    cs = c_hist_unconn[0,:]    
    plt.figure(fig_nr+2)        
    plt.clf()
    #Cs = np.linspace(-2,1)
    #plt.plot(Cs,-np.log((1-np.exp(Cs))/(2)))
    numblocks = len(c_hist_unconn[1:,0])
    Czeros = np.zeros(numblocks)    
    scale = 1e6
    for i in range(1,numblocks+1):
    #for i in range(1,2):
        b1u = c_hist_unconn[i,:]/scale
        b1c = c_hist_conn[i,:]/scale
        #diff = b1u-b1c
        diffl = np.log(b1u)-np.log(b1c)
        poly = np.polyfit(cs,diffl,1)
        C = -poly[1]/poly[0]
        if 1:        
            #plt.plot(cs,np.log(b1u),'b')
            #plt.plot(cs,np.log(b1c),'r')
            plt.plot(cs,diffl,'c')
            plt.plot(cs,np.polyval(poly,cs),'x--')
            plt.plot([C],[0],'bo')
        Czeros[i-1]=C
    #Czeros *= -1
    #print(Czeros)
    if(opt=='FB'):
        dF = -np.log((1-np.exp(Czeros))/(1+np.exp(Czeros))/deg)
        dF_a = -np.log((1-np.exp(np.mean(Czeros)))/(1+np.exp(np.mean(Czeros)))/deg)
    if(opt=='FD'):
        dF = -np.log((1-np.exp(Czeros))/2)
        dF_a = -np.log((1-np.exp(np.mean(Czeros)))/2)
    if(opt=='mass'):
        dF = -Czeros/beta
        dF_a = -np.mean(Czeros)/beta
    print("Option: "+opt)
    print("F_O-F_oo="+str(np.mean(-Czeros/beta))+"+/-"+str(np.std(Czeros/beta)))
    print("dF="+str(np.mean(dF/beta))+"+/-"+str(np.std(dF/beta)))
    print("dF_alternative_mean="+str(dF_a/beta))    
    avg_sign = (1-np.exp(Czeros))/(1+np.exp(Czeros))
    print("<sign>="+str(np.mean(avg_sign))+"+/-"+str(np.std(avg_sign))+'\n')
    
def free_energy_diff_3p(f,fig_nr,beta=1.0,P=10,opt='FB',deg=1):    
    c01 = np.loadtxt(f+'diag0-1/fsum_N2_P'+str(P)+'.dat')[:,1:]
    c10 = np.loadtxt(f+'diag1-0/fsum_N1_P'+str(P)+'.dat')[:,1:]
    c02 = np.loadtxt(f+'diag0-2/fsum_N2_P'+str(P)+'.dat')[:,1:]
    c20 = np.loadtxt(f+'diag2-0/fsum_N1_P'+str(P)+'.dat')[:,1:]
    cs = c01[0,:]
    numblocks = len(c01[1:,0])
    plt.figure(fig_nr)        
    plt.clf()
    #Cs = np.linspace(-2,1)
    #plt.plot(Cs,-np.log((1-np.exp(Cs))/(2)))
    Czeros1 = np.zeros(numblocks)    
    Czeros2 = np.zeros(numblocks)    
    scale = 1e6
    for i in range(1,numblocks+1):
        b01 = c01[i,:]/scale
        b10 = c10[i,:]/scale
        b02 = c02[i,:]/scale
        b20 = c20[i,:]/scale
        diffl1 = np.log(b01)-np.log(b10)
        diffl2 = np.log(b02)-np.log(b20)
        poly1 = np.polyfit(cs,diffl1,1)
        poly2 = np.polyfit(cs,diffl2,1)
        C1 = -poly1[1]/poly1[0]
        C2 = -poly2[1]/poly2[0]
        if 1:        
            plt.plot(cs,np.log(b01),'b')
            plt.plot(cs,np.log(b10),'r')
            plt.plot(cs,diffl1,'c')
            plt.plot(cs,np.polyval(poly1,cs),'x--')
            plt.plot([C1],[0],'bo')
        Czeros1[i-1]=C1
        Czeros2[i-1]=C2
    print("F_Oo-F_ooo="+str(np.mean(Czeros1)))
    print("F_D-F_ooo="+str(np.mean(Czeros2)))    
    if(opt=='FB'):
        dF = -np.log((1-3*np.exp(Czeros1)+2*np.exp(Czeros2))/(1+3*np.exp(Czeros1)+2*np.exp(Czeros2))/deg)/beta
        #dF = -np.log((1-np.exp(Czeros1)+np.exp(Czeros2))/(1+np.exp(Czeros1)+np.exp(Czeros2)))/beta
        dF2 = -np.log((1-3*np.exp(Czeros1))/(1+3*np.exp(Czeros1))/deg)/beta
        sign = (1-3*np.exp((Czeros1))+2*np.exp((Czeros2)))/(1+3*np.exp((Czeros1))+2*np.exp((Czeros2)))/deg
        print("sign="+str(np.mean(sign))+"+/-"+str(np.std(sign)))
        dF_a = -np.log((1-3*np.exp(np.mean(Czeros1))+2*np.exp(np.mean(Czeros2)))/(1+3*np.exp(np.mean(Czeros1))+2*np.exp(np.mean(Czeros2)))/deg)/beta
    if(opt=='FD'):
        dF = -np.log((1-3*np.exp(Czeros1)+2*np.exp(Czeros2))/6)/beta
        dF_a = -np.log((1-3*np.exp(np.mean(Czeros1))+2*np.exp(np.mean(Czeros2)))/6)/beta
    if(opt=='spin1half'):
        dF = -np.log((1-0*np.exp(Czeros1)-np.exp(Czeros2))/(1+0*np.exp(Czeros1)+np.exp(Czeros2))/deg)/beta
        dF_a = -np.log((1-0*np.exp(np.mean(Czeros1))-np.exp(np.mean(Czeros2)))/(1+0*np.exp(np.mean(Czeros1))-np.exp(np.mean(Czeros2)))/deg)/beta
        print(((1-np.exp(np.mean(Czeros2)))/(1+np.exp(np.mean(Czeros2)))/deg)/beta)
        
    print("Option: "+opt)
    print("dF="+str(np.mean(dF))+"+/-"+str(np.std(dF)/np.sqrt(numblocks)))        
    #print("dF_twodiags="+str(np.mean(dF2))+"+/-"+str(np.std(dF2)/np.sqrt(numblocks)))
    print("dF_alternative_mean="+str(dF_a) +'\n')     
    
def plot_theoretical_energies(fig_nr,clear=False,d=1,hw=3.0,color='k',label='',ax=None):
    T = np.linspace(1.5,100,200)
    kB = 1/11.6045
    beta = 1.0/(kB*T)
    tau = 0.1
    tmp = beta/tau
    #P = map(lambda i: int(i), tmp)
    P = tmp.astype(int)
    #print(P)
    for i,p in enumerate(P):
        if(np.abs(tmp[i]-p)<0.01):
            P[i] += 1
        if(p<2):
            P[i] = 2
    #print(P)
    #P = np.floor(beta/tau)
    hwb = hw*beta
    x = hwb/P
    b = 1 + 0.5*x**2 + 0.5*x*np.sqrt(4+x**2)
    bder = x + 0.5*np.sqrt(4+x**2) + 0.5*x**2/np.sqrt(4+x**2)
    frac1 = (2*np.sinh(hwb/2)**2/np.sinh(hwb))**d
    frac2 = np.tanh(hwb/2)**(-1)
    frac3 = np.tanh(hwb)**(-1)
    const = 1
    E = 2*d*(0.5+1/(np.exp(hwb)-1))
    if 1:
        frac1 = ((b**P-1)**2/(b**(2*P)-1))**d
        frac2 = (b**P+1)/(b**P-1)
        frac3 = (b**(2*P)+1)/(b**(2*P)-1)
        const = bder/b
        E = 2*const*d*(0.5+1/(b**P-1)) #distinguish
        if 0:
            for i,p in enumerate(P):
                if(p==1):
                    E[i]=2*d/hwb[i]
    Eb = const*d*(frac2 + frac3*frac1)/(1+frac1) #boson
    Ef = const*d*(frac2 - frac3*frac1)/(1-frac1) #fermion
    if 0:    
        fileb = open(files[0],'w')
        filef = open(files[1],'w')
        for i in range(len(Eb)):
            fileb.write(str(T[i])+"\t"+str(Eb[i]*hw)+"\n")
            filef.write(str(T[i])+"\t"+str(Ef[i]*hw)+"\n")
    #Results are in units of hw
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    if ax is None:
        plt.plot(T,E,'--',color=color,label=label)
        plt.plot(T,Eb,'--',color=color)
        plt.plot(T,Ef,'--',color=color)
    else:
        ax.plot(T,E,'--',color=color,label=label)
        ax.plot(T,Eb,'--',color=color)
        ax.plot(T,Ef,'--',color=color)
    
def plot_bennett_vs_T(f,fig_nr,clear=True):
    bennett = np.loadtxt(f+'energies_bennett.dat',comments='%')
    fermion = np.loadtxt(f+'energies_fermion.dat',comments='%')
    beta = bennett[:,0]
    Tb = 11.6045/beta
    en_scale = hw
    b = bennett[:,1]/en_scale
    berr = bennett[:,2]/en_scale
    DF = bennett[:,3]/en_scale
    DFerr = bennett[:,4]/en_scale
    sgn_ben = bennett[:,5]
    sgnerr_ben = bennett[:,6]
    WB = bennett[:,7]
    WBerr = bennett[:,8]
    beta2 = fermion[:,0]
    Tf = 11.6045/beta2
    f = fermion[:,1]/en_scale
    ferr = fermion[:,2]/en_scale
    WF = fermion[:,3]
    WFerr = fermion[:,4]
    sgn_expavg = WF/WB  
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    plt.plot(Tb,b+(-np.log(sgn_expavg)/beta)/en_scale,'bo',label='Exponent average')
    plt.errorbar(Tb,b+DF,berr+DFerr,color='r',linestyle='',marker='s',label='Bennett method')
    plt.errorbar(Tf,f,ferr,color='g',linestyle='',marker='v',label='Direct method')
    plt.xlabel('$T~\mathrm{(K)}$')
    plt.ylabel(r'$E/\hbar\omega_0$')
    plt.legend(loc='upper left',fontsize=22)
    plt.title('Fermions at low temperature')
    plt.xlim([0,26])
    plt.ylim([1.85,2.5])
    
    plt.figure(fig_nr+1)    
    plt.clf()
    ind = np.isfinite(beta) & np.isfinite(sgn_expavg) & (sgn_expavg>0)
    ind2 = np.isfinite(beta) & np.isfinite(sgn_ben) & (sgn_ben>0)
    poly1,cov1 = np.polyfit(beta[ind],np.log(sgn_expavg[ind]),1,cov=True)
    print(poly1)
    poly2,cov2 = np.polyfit(beta[ind2],np.log(sgn_ben[ind2]),1,cov=True)
    print(poly2)
    plt.plot(beta,sgn_expavg,'go',label='Exp. avg: 3.309 meV')
    plt.plot(beta,sgn_ben,'rs',label='Bennett: 2.976 meV')
    plt.xlabel(r'$\beta\,(\mathrm{meV}^{-1})$')
    x = np.linspace(beta[0],beta[-2])
    plt.plot(x,np.exp(np.polyval(poly1,x)),'g--')
    plt.plot(x,np.exp(np.polyval(poly2,x)),'r--')
    ax = plt.gca()
    ax.set_yscale('log')
    plt.ylabel(r'$\log\,\langle\mathrm{sgn}\rangle$')
    plt.legend(loc='upper right',fontsize=22)
    plt.xlim([0.4,1.6])
    plt.ylim([0.005,0.5])
    ax.tick_params(length=6,which='major')
    ax.tick_params(length=3,which='minor')
    plt.title('1D')

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
    f10= '../run10/'
    f11= '../run11/'
    f12= '../run12/'
    f13= '../run13/'
    f14= '../run14/'
    f15= '../run15/'
    f16= '../run16/'
    f17= '../run17/'
    f18= '../run18/'
    f19= '../run19/'
    f20= '../run20/'
    f21= '../run21/'
    f22= '../run22/'
    f23= '../run23/'
    s=['s0-25/','s0-5/','s1/','s2/','s4/']

    
    
    fi1 ='../ideal/boson/1D/beta1/MetaD/disconnected/'
    fi2 ='../coulomb_extra3/RW1-34/disconnected/'
    fi3 = '../coulomb_extra3/RW1-34/connected/'
    fi2a ='../coulomb/circular/RW1-4/singlet/beta1/disconnected/'
    fi3a = '../coulomb/circular/RW1-4/triplet/beta1/'
    fi2b ='../coulomb/circular/RW2/singlet/disconnected/'
    fi3b = '../coulomb/circular/RW2/triplet/'
    fi4 ='../ideal/fermion/1D/beta1/MetaD/'
    fiucon = '../ideal/boson/2D/beta1-5/disconnected/' #'../three/spin3half/beta0-5/MetaD_CV8/diag0-1/'
    ficonn = '../ideal/boson/2D/beta1-5/connected/'#'../three/spin3half/beta0-5/MetaD_CV8/diag0-1_reversed/'
    fjucon = '../coulomb_extra4/RW1/disconnected/'   
    fjconn = '../coulomb_extra4/RW1/connected/'   
    fi5 ='../coulomb/anisotropy1-38/singlet/beta1/MetaD/disconnected/'
    fi6 ='../coulomb/anisotropy1-38/triplet/beta1/MetaD/'
    fi5b ='../coulomb/anisotropy3/singlet/beta1/disconnected/'
    fi6b ='../coulomb/anisotropy3/triplet/beta1/'
    f_D1 ='../ideal/distinguish/1D/Energy_vs_T/'
    f_B1 ='../ideal/boson/1D/Energy_vs_T/'
    f_F1 ='../ideal/fermion/1D/Energy_vs_T/'
    f_D2 ='../ideal/distinguish/2D/Energy_vs_T/'
    f_B2 ='../ideal/boson/2D/Energy_vs_T/'
    f_F2 ='../ideal/fermion/2D/Energy_vs_T/'
    f_FwoM='../ideal/fermion/1D/beta1/woMetaD/'
    f_FwMshort='../ideal/fermion/1D/beta1/MetaD/'
    f_FwM='../ideal/fermion/1D/beta1/MetaD_longer/'
    f_BwoM='../ideal/boson/1D/beta1/woMetaD/disconnected/'
    f_BwM='../ideal/boson/1D/beta1/MetaD_longer/disconnected/'
    f_FwoM2='../ideal/fermion/2D/beta1/woMetaD/'
    f_FwM2='../ideal/fermion/2D/beta1/MetaD/'
    f_BwoM2='../ideal/boson/2D/beta1/woMetaD/'
    f_BwM2='../ideal/boson/2D/beta1/MetaD/disconnected/'
    fi8 ='../ideal/distinguish/1D/beta1/Energy_vs_P/'
    fi9 ='../ideal/distinguish/1D/beta2/Energy_vs_P/'
    fi10='../ideal/fermion/2D/beta1/woMetaD_longer/'
    fid ='../ideal/distinguish/1D/beta1/'
    fid2 ='../ideal/distinguish/2D/beta1/'
    
    
    fbennett1D = '../ideal/fermion/1D/Energy_vs_T/FreeEnergyDiff/'
    fbennett2D = '../ideal/fermion/2D/Energy_vs_T/FreeEnergyDiff/'
    fbennett2 = '../coulomb/anisotropy1-38/triplet/Energy_vs_T/FreeEnergyDiff/'
    f_S1 = '../coulomb/anisotropy1-38/singlet/beta1/MetaD/disconnected/'
    f_T1 = '../coulomb/anisotropy1-38/triplet/beta1/MetaD/'
    f_S = '../coulomb/anisotropy1-38/singlet/Energy_vs_T/disconnected/'
    f_T = '../coulomb/anisotropy1-38/triplet/Energy_vs_T/'
    f3p = '../three/spin3half/beta0-5/'


    fr2sing = '../coulomb170614/circular/RW4/singlet/'
    fr2trip = '../coulomb170614/circular/RW1-4/triplet/beta1/'
    #plot_bennett_vs_T(fbennett,1)
    #plot_theoretical_energies(1,0)
    
    #free_energy_diff(fjucon,fjconn,5,beta=0.5,P=13,opt='FB',deg=2)
        
    if 0:
        plot_rAB(fi6,2,1,2,'r',marker='x',name='Disting.',linestyle='-',show_errors=1)
        plot_rAB(fi7,2,0,2,'b',marker='x',name='Boson',linestyle='-',show_errors=1)
        plot_rAB(fi8,2,0,2,'g',marker='x',name='Fermion',linestyle='-',show_errors=1)
        plt.legend(loc='upper right',fontsize=22)    
        plt.title(r'$\sigma_\mathrm{LJ}=0.5\,\mathrm{nm}$')
        plot_gauss_data(fi7,3,'Wt')
        plot_s_int(fi8,4,1,'')
        
    if 0:
        plot_2d_dist([fi5,fi6],1,[r'$\mathrm{Singlet}$',r'$\mathrm{Triplet}$',r'$\mathrm{Partial}$'],
                     '$\mathrm{Anisotropy~}\eta = 1.38$',use_contour=True,stride=5,num_smooths=0)

    if 0:
        plot_2d_dist([fi5b,fi6b],2,[r'$\mathrm{Singlet}$',r'$\mathrm{Triplet}$',r'$\mathrm{Partial}$'],
                     '$\mathrm{Anisotropy~}\eta = 3.0$',use_contour=True,stride=5,num_smooths=0)

    if 0:
        plot_2d_dist([fi1,fi4],1,[r'$\mathrm{Boson}$',r'$\mathrm{Fermion}$',r'$\mathrm{Partial}$'],
                     r'$\mathrm{With~Metadynamics}$',use_contour=True,stride=5,num_smooths=3)
   
    if 0:
        plot_2d_dist_temps('../coulomb/anisotropy1-38/',9,0)

    if 0:
        plot_2dpaircorr([fi5b,fi6b],2,[r'$\mathrm{Singlet}$',r'$\mathrm{Triplet}$'],
                        r'$\mathrm{Pair~correlation,~}\eta = 3.0$',num_smooths=15)

    if 0:
        plot_2dpaircorr([fi5,fi6],4,[r'$\mathrm{Singlet}$',r'$\mathrm{Triplet}$'],
                        r'$\mathrm{Pair~correlation,~}\eta = 1.38$',num_smooths=15,m=2)
        plt.xlim([0,70])

        
    if 0:
        plot_rAB(fi2a,0,True,2,'b',scale=100,marker='s',name=r'$R_\mathrm{W}=1.4,~\mathrm{Singlet}$',linestyle='-',show_errors=2)
        plot_rAB(fi3a,0,False,2,'g',scale=100,marker='D',name=r'$R_\mathrm{W}=1.4,~\mathrm{Triplet}$',linestyle='-',show_errors=2)
        plot_rAB(fi2b,0,False,2,'r',scale=100,marker='^',name=r'$R_\mathrm{W}=2.0,~\mathrm{Singlet}$',linestyle='-',show_errors=2)
        plot_rAB(fi3b,0,False,2,'darkorange',scale=100,marker='o',name=r'$R_\mathrm{W}=2.0,~\mathrm{Triplet}$',linestyle='-',show_errors=2)
        plt.legend(loc='upper right',fontsize=22)
        plt.ylim([0,3.3])
        plt.xlim([0,85])
        plt.ylabel(r'$g(r)~(10^{-2}\,\mathrm{nm}^{-1})$')
        plt.xlabel(r'$r~(\mathrm{nm})$')
        plt.title(r'$\mathrm{Radial~pair~correlation,~circular}$')

                                
    if 0:
        r,_,_ = plot_rAB(fi1,0,True,2,'b',marker='o',name=r'Boson',linestyle='',show_errors=4)
        r2,_,_ = plot_rAB(fi4,0,False,2,'g',marker='D',name=r'Fermion',linestyle='',show_errors=4)
        r3,_,_ = plot_rAB(fid,0,False,2,'r',marker='s',name=r'Distinguishable',linestyle='',show_errors=8)
        plot_rAB_th(0,r,2,'bos')        
        plot_rAB_th(0,r2,2,'fer')
        plot_rAB_th(0,r3,2,'dis')        
        plt.legend(loc='upper right',fontsize=22)
        plt.title(r'$\mathrm{2D,~with~Metadynamics}$')
        #plt.ylim([-0.005,0.025])
        plt.xlim([0,20])
        
        #plot_rAB(f,fig_nr,clear=True,d=2,color='blue',marker='x',name=None,linestyle='-',show_errors=1):
    #plot_rAB(fr2sing,1,1)
    #plot_rAB(fr2trip,1,0,color='red')
        
    #free_energy_diff(fi2,fi3,7,beta=1.0,P=15,deg=2)       
    #plt.title('1D')
        
    #plot_s_int(fi2,3,1,'')
    #plot_s_int(fi3,4,1,'')
        

    
    if 0:        
        #plot_cv('../coulomb/circular/RW3/singlet/woMetaD/shortrun/',2)
        plot_s_int('../coulomb/circular/RW3/singlet/',4,1,'log',color='b',labels=['With MetaD',''],linestyle='-.')
        plot_s_int('../coulomb/circular/RW3/singlet/woMetaD/',1,1,'log',color='b',labels=['Without MetaD',''],linestyle='-.')
        #plot_s_int('../coulomb/circular/RW3/singlet/woMetaD/',4,0,'log',labels=['With wall, rew.','With wall, part. rew.'],with_wall=True,linestyle='-')
        plt.ylim([0,40])
        plt.xlim([-35,100])
        #plt.title('Free energy surface')
        plt.xlabel(r'$s=\beta\Delta U$')
        plt.ylabel(r'$F(s)~(\mathrm{meV})$')
        plot_gauss_data('../test/',6,1,opt='Ws')
        plt.title('Without wall')
        plt.xlim([-30,50])
        #plot_gauss_data('../test2/',7,1,opt='Ws')
        #plt.title('With wall')
        
    if 0:
        #plot_cv('../test5/',1,100000)  
        plot_s_int('../test3/',2,1,'')
        plt.xlim([-1,20])
        #plt.ylim([-0.5,30])
        plt.xlabel('$s$')
        plt.title('Without Metadynamics')
        plt.legend(loc='upper right',fontsize=24)
        plot_cv('../test2/',3,100000)
        plot_s_int('../test2/',4,1,'')
        plt.xlim([-1,20])
        #plt.ylim([-0.5,30])
        plt.xlabel('$s$')
        plt.title('With Metadynamics')
        plt.legend(loc='upper right',fontsize=24)
        plot_gauss_data('../test2/',5,1,opt='Wt')
        
    if 0:
        #plot_cv('../sign_cv/woMetaD/',0,100000,'')
        #plot_cv('../test3/',1,100000,'onlyperm')
        #plot_energies_vs_t(f1,0,n=100000)
        #plot_gauss_data('../test3/',2,1,'st',name='test3')
        #plot_gauss_data('../ideal/boson/2D/beta1/MetaD/disconnected/',9,1,'Ws',name='test3fast')
        """plot_s_int(fi2,2,1,'')
        #plt.ylim([-5,5])
        #plt.xlim([-1.1,1.1])
        plt.title('sing')
        plot_s_int(fi3,3,1,'')
        #plt.ylim([-5,5])
        #plt.xlim([-1.1,1.1])
        plt.title('trip')"""
        """plot_s_int('../ideal/boson/1D/beta1/MetaD/disconnected/',5,1,'')
        #plt.ylim([-1,1])
        plt.xlim([-40,40])
        plt.title('boson')"""
        """plot_s_int('../test4/',6,1,'')
        plt.xlim([-40,40])
        plt.title('gaussians for force')"""
    
    #W_distr
    if 0:
        #Ws,pW = plot_cont(f_FwoM ,4,-1,100000,cv='G')
        #Ws2,pW2 = plot_cont(f_FwMshort,5,-1,100000,cv='G')
        plt.figure(6)
        plt.clf()
        maxi=192
        plt.plot(Ws,pW,'b',label='Without Metadynamics',linewidth=4)
        plt.plot(Ws2[:maxi:5],pW2[:maxi:5],'r',label='With Metadynamics',linestyle='',marker='.',linewidth=4,markersize=13)
        plt.plot(Ws2[maxi::],pW2[maxi::],'r',linestyle='',marker='.',markersize=13)
        plt.xlim([-480,40])
        #ax = plt.gca()
        #ax.set_yscale('log')
        plt.legend(loc='upper left',fontsize=24)
        #plt.title('Fermion')
        plt.xlabel('$W$')
        plt.ylabel(r'$\log\,p(W)$')
        plt.tight_layout()
        
    
    #New 1D pair corr plot
    if 0:
        fig = plt.figure(num=3,figsize=(12,4))
        gs = mpl.gridspec.GridSpec(1,2)
        gs.update(hspace=0.0, wspace=0.1, right=0.85, bottom=0.22, top=0.95)
        plt.clf()

        ax0 = plt.subplot(gs[0])
        d=1
        fig_nr=3
        n=10
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        r,pb,pberr=plot_rAB(f_BwM,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n)
        r,_,_=plot_rAB(f_FwM,fig_nr,0,d,'g','D',name=r'Fermion',linestyle='',show_errors=n)
        #plt.text(1.5, 2.4, r'1D',fontsize=40)
        plt.text(15, 2.6, r'1D',fontsize=35)
        plt.xlim([0,20])
        plt.ylim([0,6])
        plt.xlabel('$r~(\mathrm{nm})$')
        plt.ylabel(r'$g(r)~(\mathrm{nm}^{-1})$')
        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.79,0.74,0.15,0.15])            
            ni=10
            m=32
            rins = r[1:m:ni]
            pins = pb[1:m:ni]
            axins.plot(r[1:m],pbth[1:m],'b')
            axins.errorbar(rins,pb[1:m:ni],pberr[1:m:ni],markersize=10,marker='o',color='b',linestyle='')        
            #axins.errorbar(rins,pd[1:22:ni],pderr[1:22:ni],markersize=10,marker='s',color='r',linestyle='')   
            #axins.plot(r[1:22],pdth[1:22],'r')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            axins.set_ylim([4.0,5.4])
            axins.set_xticks(np.array([0,1,2,3]))        
            axins.set_yticks(np.array([4.0,4.5,5.0]))
        d=2
        ax1 = plt.subplot(gs[1])
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        r,pb,pberr=plot_rAB(f_BwM2,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n2)
        r,_,_=plot_rAB(f_FwM2,fig_nr,0,d,'g','D',name=r'Fermion',linestyle='',show_errors=n2)
        #plt.legend(loc=[0.37,0.8],fontsize=20)
        plt.legend(loc='best',fontsize=24)
        plt.text(1.5, 4.2, r'2D',fontsize=40)
        plt.xlim([0,20])
        plt.ylim([0,0.7])
        plt.text(15, 0.3, r'2D',fontsize=35)
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position("right")
        plt.ylabel(r'$g(r)~(\mathrm{nm}^{-1})$')
        plt.xlabel('$r~(\mathrm{nm})$') 
        
        
    #1D pair corr plot    
    if 0:
        fig = plt.figure(num=2,figsize=(12,8))
        gs = mpl.gridspec.GridSpec(2,2)
        gs.update(hspace=0.0, wspace=0.02, bottom=0.15, top=0.9)
        plt.clf()

        ax00 = plt.subplot(gs[0,0])
        d=1
        fig_nr=2
        n=10
        r,pb,pberr=plot_rAB(f_BwoM,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n)
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        pdth = plot_rAB_th(fig_nr,r,d,'dis')
        r,pb,pberr=plot_rAB(f_BwoM,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n)
        r,_,_=plot_rAB(f_FwoM,fig_nr,0,d,'g','D',name=r'Fermion',linestyle='',show_errors=n)
        r,pd,pderr=plot_rAB(fid,fig_nr,0,d,'r','s',name=r'Distinguishable',linestyle='-',show_errors=n)    
        plt.text(5, 4.7, r'1D',fontsize=35)
        #plt.legend(loc='upper right',fontsize=20)
        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.37,0.74,0.15,0.15])            
            ni=10
            m=32
            rins = r[1:m:ni]
            pins = pb[1:m:ni]
            axins.plot(r[1:m],pbth[1:m],'b')
            axins.plot(r[1:m],pdth[1:m],'r')
            axins.errorbar(rins,pb[1:m:ni],pberr[1:m:ni],markersize=10,marker='o',color='b',linestyle='')        
            axins.errorbar(rins,pd[1:m:ni],pderr[1:m:ni],markersize=10,marker='s',color='r',linestyle='')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            axins.set_ylim([4.0,5.4])
            axins.set_xticks(np.array([0,1,2,3]))        
            axins.set_yticks(np.array([4.0,4.5,5.0]))

        ax01 = plt.subplot(gs[0,1],sharey=ax00)
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        r,pb,pberr=plot_rAB(f_BwM,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n)
        r,_,_=plot_rAB(f_FwM,fig_nr,0,d,'g','D',name=r'Fermion',linestyle='',show_errors=n)
        #plt.text(1.5, 2.4, r'1D',fontsize=40)
        plt.text(5, 4.7, r'1D',fontsize=35)
        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.79,0.74,0.15,0.15])            
            ni=10
            m=32
            rins = r[1:m:ni]
            pins = pb[1:m:ni]
            axins.plot(r[1:m],pbth[1:m],'b')
            axins.errorbar(rins,pb[1:m:ni],pberr[1:m:ni],markersize=10,marker='o',color='b',linestyle='')        
            #axins.errorbar(rins,pd[1:22:ni],pderr[1:22:ni],markersize=10,marker='s',color='r',linestyle='')   
            #axins.plot(r[1:22],pdth[1:22],'r')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            axins.set_ylim([4.0,5.4])
            axins.set_xticks(np.array([0,1,2,3]))        
            axins.set_yticks(np.array([4.0,4.5,5.0]))
        
        ax10 = plt.subplot(gs[1,0],sharex=ax00)
        d=2
        n2 = 4
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        pdth = plot_rAB_th(fig_nr,r,d,'dis')
        r,pb,pberr=plot_rAB(f_BwoM2,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n2)
        r,_,_=plot_rAB(f_FwoM2,fig_nr,0,d,'g','D',name=r'Fermion',linestyle='',show_errors=n2)
        r,pd,pderr=plot_rAB(fid2,fig_nr,0,d,'r','s',name=r'Distinguishable',linestyle='-',show_errors=n)    
        plt.legend(loc=[0.3,0.6],fontsize=20)
        #plt.text(1.5, 4.2, r'2D',fontsize=40)
        plt.xlim([0,20])
        plt.text(5, -0.05, r'2D',fontsize=35)

        
        ax11 = plt.subplot(gs[1,1],sharey=ax10,sharex=ax01)
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        r,pb,pberr=plot_rAB(f_BwM2,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n2)
        r,_,_=plot_rAB(f_FwM2,fig_nr,0,d,'g','D',name=r'Fermion',linestyle='',show_errors=n2)
        plt.legend(loc=[0.37,0.8],fontsize=20)
        plt.text(1.5, 4.2, r'2D',fontsize=40)
        plt.xlim([0,20])
        plt.ylim([-0.1,0.7])
        plt.text(5, -0.05, r'2D',fontsize=35)

        ax_main = fig.add_subplot(111,frameon=False)
        plt.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
        ax_main.set_ylabel(r'$\mathrm{Energy}~(\hbar\omega_0)$',fontsize=30)     
        plt.ylabel(r'$g(r)~(\mathrm{nm}^{-1})$')
        ax_main.set_xlabel('$r~(\mathrm{nm})$',fontsize=30)
    
        plt.setp(ax00.get_xticklabels(),visible=False)
        yticks = ax10.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.setp(ax01.get_xticklabels(),visible=False)
        yticks = ax11.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.setp(ax01.get_yticklabels(),visible=False)
        plt.setp(ax11.get_yticklabels(),visible=False)
        xticks = ax10.xaxis.get_major_ticks()
        xticks[-1].label1.set_visible(False)        
        ax00.set_title('Without Metadynamics',fontsize=30)
        ax01.set_title('With Metadynamics',fontsize=30)
        
        
 

        #plt.ylim([-0.02,0.1])
        #plt.ylim([0,0.2])
        
            #plt.xticks(visible=False)
        #plt.yticks(visible=False)
        #mark_inset(ax,axins,loc1=1,loc2=4,fc="none",ec="0.5")
        #fig.sca(ax)       
       
        
        
        
    #1D density plot
    if 0: 
        """fig = plt.figure(num=0,figsize=(12,5))
        gs = mpl.gridspec.GridSpec(1,2)
        gs.update(wspace=0.05,bottom=0.15,top=0.9,left=0.16,right=0.95)
        plt.clf()
        r=plot_1d_dist([f_BwoM,f_FwoM],0,0,n=5,ax=ax0)
        plt.clf()
        ax0 = plt.subplot(gs[0])
        plot_1d_dist_th(0,r,d=1,hw=3.0,ax=ax0)        
        r=plot_1d_dist([f_BwoM,f_FwoM],0,0,n=5,ax=ax0)
        plt.title('Without Metadynamics',fontsize=30)
        plt.xlabel(r'$r~(\mathrm{nm})$')
        plt.ylabel(r'$p(r)$')"""
    
        fig = plt.figure(num=0,figsize=(10,6))
        plt.clf()
        ax1 = plt.gca()
        #ax1 = plt.subplot(gs[1],sharey=ax0)
        plot_1d_dist_th(0,r,d=1,hw=3.0,ax=ax1)
        r = plot_1d_dist([f_BwM,f_FwM],0,0,n=5,ax=ax1)
        #plt.title('Probability density',fontsize=30)
        #plt.setp(ax1.get_yticklabels(),visible=False)
        #xticks = ax0.xaxis.get_major_ticks()
        #xticks[-1].label1.set_visible(False)
        plt.xlabel(r'$r~(\mathrm{nm})$')
        #plt.legend(loc=(-0.3,0.5),fontsize=22)
        plt.legend(loc=[0.71,0.66],fontsize=22)
        plt.ylabel(r'$p(r)$')
        
    #energies and paircorr
    if 0:
        fig = plt.figure(num=1,figsize=(15,8))   
        gs = mpl.gridspec.GridSpec(2,2,height_ratios=[1,1])
        gs.update(hspace=0.0,wspace=0.07, bottom=0.12,top=0.95,left=0.06,right=0.9)
        plt.clf()
        ax0 = plt.subplot(gs[0])
        plot_theoretical_energies(1,False,1,3.0,'k',ax=ax0)#,r'$\mathrm{Theory}$')
        plot_energies(f_B1,1,0,'beta',ax=ax0,label=r'Boson',marker='D',color='b')
        plot_energies(f_F1,1,0,'beta',ax=ax0,label=r'Fermion',marker='s',color='g')
        plot_energies(f_D1,1,0,'beta',ax=ax0,label=r'Distinguishable',marker='o',color='r')
        plt.xlim([0,80])
        plt.ylim([0,6])
        plt.ylabel(r'$E/\hbar\omega_0$')
        plt.text(4, 4.8, r'1D',fontsize=35)
        #plt.plot([0,10],[2,2],'k--',label=r'$\mathrm{Theoretical~value}$')
        ax1 = fig.add_subplot(gs[2],sharex=ax0)
        plot_theoretical_energies(1,False,2,3.0,'k',ax=ax1)#,r'$\mathrm{Theory}$')
        plot_energies(f_B2,1,0,'beta',ax=ax1,label=r'Boson',marker='D',color='b')
        plot_energies(f_F2,1,0,'beta',ax=ax1,label=r'Fermion',marker='s',color='g')
        plot_energies(f_D2,1,0,'beta',ax=ax1,label=r'Distinguishable',marker='o',color='r')
        plt.xlim([0,80])
        plt.ylim([0,10.])
        plt.setp(ax0.get_xticklabels(),visible=False)
        yticks = ax1.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        #plt.legend(loc=[0.2,0.73],fontsize=23)
        plt.text(4, 8, r'2D',fontsize=35)
        plt.xlabel(r'$T\,(\mathrm{K})$')
        plt.ylabel(r'$E/\hbar\omega_0$')    
        
        ax2 = fig.add_subplot(gs[1])
        d=1
        fig_nr=1
        n=10
        r,pb,pberr=plot_rAB(f_BwoM,fig_nr,0,d,'b','o',name=r'Boson',linestyle='',show_errors=n)
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        pdth = plot_rAB_th(fig_nr,r,d,'dis')
        r,pb,pberr=plot_rAB(f_BwoM,fig_nr,0,d,'b','D',name=r'Boson',linestyle='',show_errors=n)
        r,_,_=plot_rAB(f_FwoM,fig_nr,0,d,'g','s',name=r'Fermion',linestyle='',show_errors=n)
        r,pd,pderr=plot_rAB(fid,fig_nr,0,d,'r','o',name=r'Distinguishable',linestyle='',show_errors=n)    
        plt.text(4, 4.7, r'1D',fontsize=35)
        plt.ylim([0,6])
        #plt.legend(loc='upper right',fontsize=20)
        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.74,0.73,0.15,0.2])            
            ni=10
            m=32
            rins = r[1:m:ni]
            pins = pb[1:m:ni]
            axins.plot(r[1:m],pbth[1:m],'b')
            axins.plot(r[1:m],pdth[1:m],'r')
            axins.errorbar(rins,pb[1:m:ni],pberr[1:m:ni],markersize=10,marker='D',color='b',linestyle='')        
            axins.errorbar(rins,pd[1:m:ni],pderr[1:m:ni],markersize=10,marker='o',color='r',linestyle='')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            axins.set_ylim([4.0,5.4])
            axins.set_xticks(np.array([0,1,2,3]))        
            axins.set_yticks(np.array([4.0,4.5,5.0]))

        ax3 = plt.subplot(gs[3],sharex=ax2)
        d=2
        n2 = 4
        pbth = plot_rAB_th(fig_nr,r,d,'bos')
        plot_rAB_th(fig_nr,r,d,'fer')
        pdth = plot_rAB_th(fig_nr,r,d,'dis')
        r,pb,pberr=plot_rAB(f_BwoM2,fig_nr,0,d,'b','D',name=r'Boson',linestyle='',show_errors=n2)
        r,_,_=plot_rAB(f_FwoM2,fig_nr,0,d,'g','s',name=r'Fermion',linestyle='',show_errors=n2)
        r,pd,pderr=plot_rAB(fid2,fig_nr,0,d,'r','o',name=r'Distinguishable',linestyle='',show_errors=n)    
        plt.legend(loc=[0.36,0.5],fontsize=24)
        #plt.text(1.5, 4.2, r'2D',fontsize=40)
        plt.xlim([0,20])
        plt.ylim([0,0.7])
        plt.text(4, 0.56, r'2D',fontsize=35)
        
        plt.setp(ax2.get_xticklabels(),visible=False)
        yticks = ax3.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        ax2.yaxis.tick_right()
        ax3.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax3.yaxis.set_label_position("right")
        ax2.set_ylabel(r'$g(r)~(\mathrm{nm}^{-1})$')
        ax3.set_ylabel(r'$g(r)~(\mathrm{nm}^{-1})$')
        ax3.set_xlabel('$r~(\mathrm{nm})$',fontsize=30)
        
        
    # FES and Wps
    if 0:
        fig = plt.figure(num=1,figsize=(14,8))   
        gs = mpl.gridspec.GridSpec(2,2,width_ratios=[1,1.5])
        gs.update(hspace=0.0,wspace=0.07, bottom=0.12,top=0.95,left=0.1,right=0.89)
        plt.clf()
        ax0 = plt.subplot(gs[:,0])  
        dE_unconn = np.loadtxt('../ideal/boson/1D/beta1/woMetaD_longer/disconnected/CV_distributions.dat')[:,0:2]
        dE_conn = np.loadtxt('../ideal/boson/1D/beta1/woMetaD_longer/connected/CV_distributions.dat')[:,0:2]
        #print(c_hist_conn)
        s = dE_unconn[:,0]
        ds = s[1]-s[0]
        hist_oo = dE_unconn[:,1]
        hist_oo /= sum(hist_oo)*ds
        hist_O = dE_conn[:,1]
        hist_O /= sum(hist_O)*ds
        maxi=911
        ax0.plot(s[maxi:], hist_oo[maxi:], 'b', linewidth=4, label='Without MetaD')
        dE_unconn = np.loadtxt('../ideal/boson/1D/beta1/MetaD_longer/disconnected/CV_distributions.dat')[:,0:2]
        hist_oo = dE_unconn[:,1]
        hist_oo /= sum(hist_oo)*ds
        ax0.plot(s[:1000:5], hist_oo[:1000:5],'r.',markersize=13,linestyle='',label='With MetaD')
        ax0.plot(s[1000::20], hist_oo[1000::20],'r.',markersize=13,linestyle='',label='With MetaD')
        ax0.set_yscale('log')
        plt.ylim([1e-9,1])
        plt.xlim([-20,20])
        plt.ylabel(r'$p(s)$')
        plt.xlabel(r'$s=\beta(U_O-U_{oo})$')
        plt.text(5, 2e-9, r'(a)',fontsize=40)
        #ax0.yaxis.tick_right()
        #ax0.yaxis.set_label_position("right")
        """ax1 = plt.subplot(gs[2],sharex=ax0)  
        dE_unconn = np.loadtxt('../ideal/boson/1D/beta1/MetaD_longer/disconnected/CV_distributions.dat')[:,0:2]
        dE_conn = np.loadtxt('../ideal/boson/1D/beta1/MetaD_longer/connected/CV_distributions.dat')[:,0:2]
        #print(c_hist_conn)
        s = dE_unconn[:,0]
        ds = s[1]-s[0]
        hist_oo = dE_unconn[:,1]
        hist_oo /= sum(hist_oo)*ds
        hist_O = dE_conn[:,1]
        hist_O /= sum(hist_O)*ds
        ax1.plot(s, -np.log(hist_oo), 'g', label='Disc.')
        ax1.plot(s, -np.log(hist_O),'k--',label='Conn.' )
        plt.setp(ax0.get_xticklabels(),visible=False)
        yticks = ax1.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)        
        plt.xlim([-60,60])
        plt.ylim([0,20])
        plt.xlabel(r'$s$')
        plt.ylabel(r'$F(s) = - \log\,p(s)$')
        plt.text(-55, 1, r'With MetaD',fontsize=30)
        plt.legend(loc=[0.65,0.85],fontsize=22)"""
            
        ax2 = plt.subplot(gs[1])
        ax2.yaxis.tick_right()
        s1,p1 = plot_s_int('../ideal/boson/1D/beta1/woMetaD_longer/disconnected/',1,0,'Whist',color='b',label=r'Without MetaD',linewidth=4,ax=ax2,withMetaD=False)
        s2,p2 = plot_s_int('../ideal/boson/1D/beta1/MetaD_longer/disconnected/',1,0,'Whist',color='r',marker='.',markersize=13,linestyle='',label=r'With MetaD',ax=ax2,withMetaD=True)
        #plt.ylim([-1,1])
        plt.xlim([-30,30])
        plt.ylim([0,1.1])
        #plt.ylabel(r'$(1-\mathrm{e}^{-s}) p(s)$')        
        plt.text(5,0.85,'Boson',fontsize=32)
        plt.text(10, 0.1, r'(b)',fontsize=40)
        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.485,0.67,0.15,0.24])            
            sins = s1[870:960]
            pins = p1[870:960]
            pins2 = p2[870:960]
            axins.plot(sins,pins,color='b',linestyle='-',linewidth=4)        
            axins.plot(sins[::7],pins2[::7],color='r',linestyle='',marker='.',markersize=13)
            axins.set_ylim([-0.005,0.05])
            #axins.errorbar(rins,pd[1:30:5],pderr[1:30:5],marker='s',color='r',linestyle='')
            #axins.plot(r[1:22],pbth[1:22],'b')
            #axins.plot(r[1:30],pdth[1:30],'r')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            #x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            plt.xticks(np.array([-12,-8,-4]))        
            plt.yticks(np.array([0,0.02,0.04]))
        ax3 = plt.subplot(gs[3],sharex=ax2)
        ax3.yaxis.tick_right()
        s1,p1 = plot_s_int('../ideal/fermion/1D/beta1/woMetaD_longer/',1,0,'Whist',color='b',ax=ax3,label='Without MetaD',linewidth=4,withMetaD=False)
        s2,p2 = plot_s_int('../ideal/fermion/1D/beta1/MetaD_longer/',1,0,'Whist',color='r',linestyle='',ax=ax3,label='With MetaD',marker='.',markersize=13,withMetaD=True)
        plt.xlim([-20,20])
        plt.ylim([-0.11,0.125])
        #ax3.legend(loc=[0.01,0.6],fontsize=22)
        ax3.legend(loc=[-0.15,0.6],fontsize=24)
        plt.setp(ax2.get_xticklabels(),visible=False)
        yticks = ax3.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.text(5,0.07,'Fermion',fontsize=32)
        plt.xlabel(r'$s=\beta(U_O-U_{oo})$')
        #ax_main = fig.add_subplot(122,frameon=False)
        #plt.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
        ax2.yaxis.set_label_position("right")
        ax3.yaxis.set_label_position("right")
        ax2.set_ylabel(r'$(1+\mathrm{e}^{-s}) p(s)$')
        ax3.set_ylabel(r'$(1-\mathrm{e}^{-s}) p(s)$')
        #ax_main.set_ylabel(r'$(1-\mathrm{e}^{-s}) p(s)$',fontsize=30)     
        #plt.tight_layout()        
        plt.text(10, -0.09, r'(c)',fontsize=40)

    #p(s)
    if 0:
        fig = plt.figure(num=4,figsize=(10,8))   
        gs = mpl.gridspec.GridSpec(2,1,height_ratios=[1,1])
        gs.update(hspace=0.0,bottom=0.12,top=0.95,left=0.16,right=0.95)
        plt.clf()
        ax0 = plt.subplot(gs[0])
        s1,p1 = plot_s_int('../ideal/boson/1D/beta1/woMetaD_longer/disconnected/',4,0,'Whist',color='b',label=r'Without MetaD',ax=ax0)
        s2,p2 = plot_s_int('../ideal/boson/1D/beta1/MetaD_longer/disconnected/',4,0,'Whist',color='r',linestyle='--',linewidth=2,label=r'With MetaD',ax=ax0)
        #plt.ylim([-1,1])
        plt.legend(loc=[0.525,0.35],fontsize=25)
        plt.xlim([-30,30])
        plt.ylim([0,1.1])
        #plt.ylabel(r'$(1-\mathrm{e}^{-s}) p(s)$')        
        plt.text(15,0.85,'Boson',fontsize=32)

        if 1:        
            #axins = fig.add_axes([0.33,0.61,0.24,0.24])
            axins = fig.add_axes([0.25,0.65,0.24,0.24])            
            sins = s1[870:960]
            pins = p1[870:960]
            pins2 = p2[870:960]
            axins.plot(sins,pins,color='b',linestyle='-')        
            axins.plot(sins,pins2,color='r',linestyle='--',linewidth=2)
            axins.set_ylim([-0.005,0.05])
            #axins.errorbar(rins,pd[1:30:5],pderr[1:30:5],marker='s',color='r',linestyle='')
            #axins.plot(r[1:22],pbth[1:22],'b')
            #axins.plot(r[1:30],pdth[1:30],'r')
            #axins.imshow(inset, extent=extent, interpolation="nearest",origin="lower")
            #x1,x2,y1,y2 = 3,10,0.15,0.25
            #axins.set_xlim(x1,x2)
            #axins.set_ylim(y1,y2)
            plt.xticks(np.array([-12,-8,-4]))        
            plt.yticks(np.array([0,0.02,0.04]))
        ax1 = plt.subplot(gs[1],sharex=ax0)
        s1,p1 = plot_s_int('../ideal/fermion/1D/beta1/woMetaD_longer/',4,0,'Whist',color='b',ax=ax1)
        s2,p2 = plot_s_int('../ideal/fermion/1D/beta1/MetaD_longer/',4,0,'Whist',color='r',linestyle='--',linewidth=2,ax=ax1)
        plt.xlim([-30,30])
        plt.ylim([-0.11,0.125])
        plt.setp(ax0.get_xticklabels(),visible=False)
        yticks = ax1.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.text(15,0.07,'Fermion',fontsize=32)
        plt.xlabel(r'$s=\beta\Delta U$')
        ax0.set_ylabel(r'$(1+\mathrm{e}^{-s}) p(s)$')
        ax1.set_ylabel(r'$(1-\mathrm{e}^{-s}) p(s)$')
        #ax_main = fig.add_subplot(111,frameon=False)
        #plt.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
        #ax_main.set_ylabel(r'$(1-\mathrm{e}^{-s}) p(s)$',fontsize=30)     
        #plt.tight_layout()
            
            



    #Ideal energies vs T
    if 0:
        fig = plt.figure(num=1,figsize=(10,8))   
        gs = mpl.gridspec.GridSpec(2,1,height_ratios=[1,1])
        gs.update(hspace=0.0,bottom=0.12,top=0.95)
        plt.clf()
        ax0 = plt.subplot(gs[0])
        plot_theoretical_energies(1,False,1,3.0,'k',ax=ax0)#,r'$\mathrm{Theory}$')
        plot_energies(f_B1,1,0,'beta',ax=ax0,label=r'Boson',marker='D',color='b')
        plot_energies(f_F1,1,0,'beta',ax=ax0,label=r'Fermion',marker='s',color='g')
        plot_energies(f_D1,1,0,'beta',ax=ax0,label=r'Distinguishable',marker='o',color='r')
        plt.xlim([0,80])
        plt.ylim([0,6])
        plt.ylabel(r'$E/\hbar\omega_0$')
        plt.text(70, 0.65, r'1D',fontsize=40)
        #plt.plot([0,10],[2,2],'k--',label=r'$\mathrm{Theoretical~value}$')
        plt.legend(loc='upper left',fontsize=22)
        ax1 = fig.add_subplot(gs[1],sharex=ax0)
        plot_theoretical_energies(1,False,2,3.0,'k',ax=ax1)#,r'$\mathrm{Theory}$')
        plot_energies(f_B2,1,0,'beta',ax=ax1,label=r'$\mathrm{Boson}$',marker='D',color='b')
        plot_energies(f_F2,1,0,'beta',ax=ax1,label=r'$\mathrm{Fermion}$',marker='s',color='g')
        plot_energies(f_D2,1,0,'beta',ax=ax1,label=r'$\mathrm{Distinguishable}$',marker='o',color='r')
        plt.xlim([0,80])
        plt.ylim([0,10.])
        plt.setp(ax0.get_xticklabels(),visible=False)
        yticks = ax1.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.text(70, 0.85, r'2D',fontsize=40)
        plt.xlabel(r'$T\,(\mathrm{K})$')
        plt.ylabel(r'$E/\hbar\omega_0$')
        
    if 0:
        fig = plt.figure(num=2,figsize=(12,8))
        gs = mpl.gridspec.GridSpec(2,2)
        gs.update(hspace=0.0, wspace=0.02, bottom=0.15, top=0.9)
        plt.clf()
        ax00 = plt.subplot(gs[0,0])
        plot_theoretical_energies(2,False,1,3.0,'k',ax=ax00)#,r'$\mathrm{Theory}$')
        plot_energies(f_B1,2,0,'beta',ax=ax00,label=r'Boson',marker='D',color='b')
        plot_energies(f_F1,2,0,'beta',ax=ax00,label=r'Fermion',marker='s',color='g')
        plot_energies(f_D1,2,0,'beta',ax=ax00,label=r'Distinguishable',marker='o',color='r')
        plt.xlim([0,25])
        plt.ylim([0.0,3])
        plt.text(1.5, 2.4, r'1D',fontsize=40)

        ax01 = plt.subplot(gs[0,1],sharey=ax00)
        plot_theoretical_energies(2,False,1,3.0,'k',ax=ax01)#,r'$\mathrm{Theory}$')        
        plot_energies(fbennett1D,2,0,'beta',ax=ax01,col=1,ediff=False,label=r'Boson',marker='D')        
        plot_energies(fbennett1D,2,0,'beta',ax=ax01,col=1,bennett=False,label=r'Fermion,direct',marker='s',color='green')
        plot_energies(fbennett1D,2,0,'beta',ax=ax01,col=1,ediff=True,label=r'Fermion, BAR',marker='o',color='darkorange')
        plt.xlim([0,25])
        plt.ylim([0.,3])
        plt.text(1.5, 2.4, r'1D',fontsize=40)

        ax10 = plt.subplot(gs[1,0])
        plot_theoretical_energies(2,False,2,3.0,'k',ax=ax10)#,r'$\mathrm{Theory}$')
        plot_energies(f_B2,2,0,'beta',ax=ax10,label=r'Boson',marker='D',color='b')
        plot_energies(f_F2,2,0,'beta',ax=ax10,label=r'Fermion',marker='s',color='g')
        plot_energies(f_D2,2,0,'beta',ax=ax10,label=r'Distinguishable',marker='o',color='r')
        plt.xlim([0,25])
        plt.ylim([1.5,5])
        plt.legend(loc=[0.3,0.8],fontsize=20)
        plt.text(1.5, 4.2, r'2D',fontsize=40)
        
        ax11 = plt.subplot(gs[1,1])
        plot_theoretical_energies(2,False,2,3.0,'k',ax=ax11)#,r'$\mathrm{Theory}$')        
        plot_energies(fbennett2D,2,0,'beta',ax=ax11,col=1,ediff=False,label=r'Boson',marker='D')        
        plot_energies(fbennett2D,2,0,'beta',ax=ax11,col=1,bennett=False,label=r'Fermion,direct',marker='s',color='green')
        plot_energies(fbennett2D,2,0,'beta',ax=ax11,col=1,ediff=True,label=r'Fermion, BAR',marker='o',color='darkorange')
        plt.xlim([0,25])
        plt.ylim([1.5,5])
        plt.legend(loc=[0.37,0.8],fontsize=20)
        plt.text(1.5, 4.2, r'2D',fontsize=40)

        ax_main = fig.add_subplot(111,frameon=False)
        plt.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
        ax_main.set_ylabel(r'$\mathrm{Energy}~(\hbar\omega_0)$',fontsize=30)     
        plt.ylabel(r'$\mathrm{Energy}~(\hbar\omega_0)$')
        ax_main.set_xlabel(r'$T\,(\mathrm{K})$',fontsize=30)
        
        plt.setp(ax00.get_xticklabels(),visible=False)
        yticks = ax10.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.setp(ax01.get_xticklabels(),visible=False)
        yticks = ax11.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)
        plt.setp(ax01.get_yticklabels(),visible=False)
        plt.setp(ax11.get_yticklabels(),visible=False)
        xticks = ax10.xaxis.get_major_ticks()
        xticks[-1].label1.set_visible(False)        
        ax00.set_title('Without Metadynamics',fontsize=30)
        ax01.set_title('With Metadynamics',fontsize=30)

    if 0:
        fig = plt.figure(num=3,figsize=(12,5))
        gs = mpl.gridspec.GridSpec(1,2)
        gs.update(wspace=0.05, bottom=0.2, top=0.95)
        plt.clf()
        bennett = np.loadtxt(fbennett1D+'energies_bennett.dat',comments='%')
        fermion = np.loadtxt(fbennett1D+'energies_fermion.dat',comments='%')
        beta = bennett[:,0]
        sgn_ben = bennett[:,5]
        sgnerr_ben = bennett[:,6]
        WB = bennett[:,7]
        WBerr = bennett[:,8]
        WF = fermion[:,3]
        WFerr = fermion[:,4]
        sgn_expavg = WF/WB  
        sgnerr_expavg = (WFerr/WF + WBerr/WB)*sgn_expavg

        ax0 = plt.subplot(gs[0])
        T = 11.6045/beta
        T[0:2]='NaN'
        plt.errorbar(T,(WF),WFerr,marker='s',color='g',markersize=12,label=r'$\langle W_\mathrm{F}\rangle$')
        plt.errorbar(T,(sgn_expavg),sgnerr_expavg,marker='D',color='b',markersize=13,label=r'$\langle \mathrm{sign}\rangle$ (exp.avg.)')
        plt.errorbar(T,(sgn_ben),sgnerr_ben,marker='o',color='r',markersize=12,label=r'$\langle \mathrm{sign}\rangle$ (BAR)')
        plt.xlabel(r'$T~(\mathrm{K})$')
        xticks = ax0.xaxis.get_major_ticks()
        xticks[-1].label1.set_visible(False)        

        ax1 = plt.subplot(gs[1])
        ax1.yaxis.tick_right()
        ax1.yaxis.set_ticks_position('both')
        beta[-1] = 'NaN'
        plt.errorbar(beta*beta[-1],(WF),WFerr,marker='s',color='g',markersize=12,label=r'$\langle W_\mathrm{F}\rangle$')
        plt.errorbar(beta,(sgn_expavg),sgnerr_expavg,marker='D',color='b',markersize=13,label=r'$\langle \mathrm{sign}\rangle$ (EA)',linestyle='')
        plt.errorbar(beta,(sgn_ben),sgnerr_ben,marker='o',color='r',markersize=12,label=r'$\langle \mathrm{sign}\rangle$ (BAR)',linestyle='')
        plt.xlabel(r'$\beta~(\mathrm{meV}^{-1})$')
        ax1.set_yscale('log')
        plt.ylim([0.003,1])
        plt.hold(True)
        legend1 = plt.legend(loc=[-0.28,0.05],fontsize=24)
        
        ind = np.isfinite(beta) & np.isfinite(sgn_expavg) & (sgn_expavg>0)
        ind2 = np.isfinite(beta) & np.isfinite(sgn_ben) & (sgn_ben>0)
        poly1,cov1 = np.polyfit(beta[ind],np.log(sgn_expavg[ind]),1,cov=True)
        print(poly1)
        print(np.sqrt(cov1[0][0]))
        poly2,cov2 = np.polyfit(beta[ind2],np.log(sgn_ben[ind2]),1,cov=True)
        print(poly2)
        print(np.sqrt(cov2[0][0]))
        #plt.plot(beta,sgn_expavg,'go',label='Exp. avg: 3.309 meV')
        #plt.plot(beta,sgn_ben,'rs',label='Bennett: 2.976 meV')
        #plt.xlabel(r'$\beta\,(\mathrm{meV}^{-1})$')
        x1 = np.linspace(beta[0]-0.1,beta[-2]+0.1)
        x2 = np.linspace(beta[2]-0.1,beta[-2]+0.1)
        l1, = plt.plot(x1,np.exp(np.polyval(poly1,x1)),'b',linestyle=':',linewidth=3)
        l2, = plt.plot(x2,np.exp(np.polyval(poly2,x2)),'r--',linewidth=2)
        #ax.tick_params(length=6,which='major')
        #ax.tick_params(length=3,which='minor')
        #plt.title('1D')

        plot_lines = [l1,l2]
        plt.legend(plot_lines, ['EA', 'BAR'], loc='best',fontsize=24)
        #plt.legend([l for l in plot_lines], ['1','2'], loc=4)
        plt.gca().add_artist(legend1)

    if 0:
        fig = plt.figure(num=3,figsize=(12,6))
        gs = mpl.gridspec.GridSpec(1,1)
        gs.update(left=0.1,bottom=0.16, top=0.95)
        plt.clf()
        ax0 = fig.add_subplot(gs[0])
        #ax0 = plt.subplot(gs[0])
        plot_energies(fbennett2,3,0,'beta',ax=ax0,col=1,ediff=False,label=r'$\mathrm{Singlet}$',marker='D',linestyle='-')        
        plot_energies(fbennett2,3,0,'beta',ax=ax0,col=1,bennett=False,label=r'$\mathrm{Triplet,direct~method}$',marker='s',color='green',linestyle='-')
        #plot_energies(f_S,3,0,'beta',ax=ax0,label=r'Singlet',marker='D',color='b',linestyle='-')
        #plot_energies(f_T,3,0,'beta',ax=ax0,label=r'Triplet',marker='s',color='g',linestyle='-')
        plt.xlim([0,90])
        plt.ylim([15,45])
        plt.xlabel(r'$T~(\mathrm{K})$')
        plt.ylabel(r'$\mathrm{Energy}~(\mathrm{meV})$')
        plt.plot([0.7,0.7],[16.35,16.35],'k<',label=r'$\mathrm{Literature~value}$',markersize=10)
        plt.plot([0.7,0.7],[18.05,18.05],'k<',markersize=10)
        axins = fig.add_axes([0.15,0.475,0.35,0.45]) 
        plot_energies(fbennett2,3,0,'beta',ax=axins,col=1,ediff=False,label=r'$\mathrm{Singlet}$',marker='D',linestyle='-')        
        plot_energies(fbennett2,3,0,'beta',ax=axins,col=1,bennett=False,label=r'$\mathrm{Triplet,direct}$',marker='s',color='green',linestyle='-')
        plot_energies(fbennett2,3,0,'beta',col=1,ediff=True,label=r"$\mathrm{Triplet,BAR}$",marker='o',color='darkorange',linestyle='-')
        axins.plot([5.2,5.2],[16.35,16.35],'k<',label=r'$\mathrm{Literature~value}$',markersize=10)
        axins.plot([5.2,5.2],[18.05,18.05],'k<',markersize=10)
        #plt.legend(loc=[1.05,0.34],fontsize=22)        
        plt.legend(loc=[1.53,-0.65],fontsize=22)        
        axins.set_xlim([5,18])        
        axins.set_ylim([15.6,20.1])        
        
        
    if 0:
        plot_energies(fbennett2D,1,1,'beta',col=1,ediff=False,label=r'$\mathrm{Singlet}$',marker='D')        
        plot_energies(fbennett2D,1,0,'beta',col=1,bennett=False,label=r'$\mathrm{Triplet,direct~method}$',marker='s',color='green')
        #plot_energies(fbennett2,1,0,'beta',col=1,ediff=True,label=r"$\mathrm{Triplet,Bennett's~method}$",marker='o',color='darkorange')
        #plt.legend(loc='upper left',fontsize=22)        
        plt.plot([0.5,0.5],[16.35,16.35],'ko',label=r'$\mathrm{Literature~value}$')
        plt.plot([0.5,0.5],[18.05,18.05],'ko')
        plt.legend(loc=[0.,0.67],fontsize=24)
        #plt.ylim([15.5,21.5])
        #plt.xlim([5,18])

    if 0:
        files = ['../ideal/boson/2D/Energy_vs_T/theory.dat','../ideal/fermion/2D/Energy_vs_T/theory.dat']
        plot_theoretical_energies(1,False,2,3.0,'k')#,r'$\mathrm{Theory}$')
        #plot_theoretical_energies(1,False,2,3.0,'k',r'$\hbar\omega_0=15\,\mathrm{meV}$')    
        #plot_theoretical_energies(1,False,2,3.0,'k',r'$\hbar\omega_0=30\,\mathrm{meV}$')    
        #plt.legend(loc='upper left',fontsize=22)
        plt.xlim([0,80])
        plt.ylim([0,6])
        plt.xlabel(r'$T\,(\mathrm{K})$')
        plt.ylabel(r'$E/\hbar\omega_0$')
        #plt.title(r'$\mathrm{2D,~without~Metadynamics}$')

    if 0:
        plot_energies(fi8,2,1,'tau',beta=1,label=r'$\beta=1~\mathrm{meV}^{-1}$',linestyle='-',marker='v')
        plot_energies(fi9,2,0,'tau',beta=2,label=r'$\beta=2~\mathrm{meV}^{-1}$',linestyle='-',marker='o',color='g')
        plt.plot([0,20],[1,1],'k--')#label=r'$\mathrm{Theoretical~value}$')
        plt.legend(loc='lower right',fontsize=26)
        
    if 0:
        file = '../coulomb/circular/wigner_scan.dat'
        file2 = '../coulomb/circular/wigner_scan_DeltaE.dat'
        data = np.loadtxt(file,comments='%')
        data2 = np.loadtxt(file2,comments='%')
        RW = data[:,0]
        plt.figure(6)
        plt.clf()
        lit = data[:,1]
        lit_err = data[:,2]
        plt.errorbar(RW,lit,yerr=lit_err,color='b',marker='D',linestyle='-',label=r'$E_\mathrm{S}~(\mathrm{lit.})$',markersize=12)
        plt.errorbar(RW,data[:,3],yerr=data[:,4],color='r',marker='o',linestyle='-',label=r'$E_\mathrm{S}~(\mathrm{calc.})$',markersize=12)
        plt.errorbar(data2[:,0],data2[:,1],yerr=(data2[:,2]),color='g',marker='s',markersize=12,linestyle='-',label=r'$E_\mathrm{ST}~(\mathrm{lit.})$')
        plt.errorbar(data2[:,0],data2[:,3],yerr=(data2[:,4]),color='darkorange',marker='^',markersize=12,linestyle='-',label=r'$E_\mathrm{ST}~(\mathrm{calc.})$')
        plt.legend(loc=[0.56,0.2],fontsize=22)
        plt.xlabel(r'$R_\mathrm{W}$')
        plt.ylabel(r'$\mathrm{Energy}~(\hbar\omega_0)$')
        plt.xlim([-0.1,2.1])
        #plt.ylim([1.5,4.5])
        #x = np.linspace(0,1.0)
        #plt.plot(x,2+x*np.sqrt(np.pi/2.0),'--')
        
    #free_energy_diff('../test/','../test3/',4,beta=1.0,P=10,opt='FB',deg=1,C=-0.101)
    #free_energy_diff('../test2/','../test4/',7,beta=1.0,P=10,opt='FB',deg=1,C=-0.488)
    free_energy_diff('../masses/oo_m1/','../masses/O_m1/',1,beta=1.0,P=10,opt='FB',deg=1,C=-0.101)
    free_energy_diff('../masses/oo_m2/','../masses/O_m2/',2,beta=1.0,P=10,opt='FB',deg=1,C=-0.5)
    free_energy_diff_mass('../masses/oo_m1/','../masses/oo_m2/','../masses/O_m1/','../masses/O_m2/',5,beta=1)  
    #print('\n')      
    #free_energy_diff_3p('../three_remote_done/beta1/',2,beta=1.0,P=10,opt='FB',deg=1)    
    
    """cv = np.loadtxt('../test6/cv.dat')
    #cv = np.loadtxt('../three/spin3half/beta0-5/woMetaD_CV5/cv.dat')
    plt.figure(0)
    plt.clf()
    plt.plot(cv[:,0],cv[:,1])
    print(np.std(cv[:,1]))"""
    
    if 1:
        plot_energies('../three/distinguish/Energy_vs_T/',1,1,'beta',color='r',linestyle='-',label='Dist')
        plot_energies('../three/boson/Energy_vs_T/',1,0,'beta',color='b',linestyle='-',marker='o',label='Boson')
        plot_energies('../three/spin3half/Energy_vs_T/',1,0,'beta',color='g',linestyle='-',marker='s',label='Fermion')
        plt.title('Three particles, 1D')
        plt.legend(loc='upper left',fontsize=22)
        plt.plot([0,10],[1.5,1.5],'k--')

    for i in plt.get_fignums():
        fig = plt.figure(i)
        fig.tight_layout()
    plt.show() 