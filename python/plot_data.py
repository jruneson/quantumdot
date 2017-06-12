import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
from matplotlib import rcParams, interactive, cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import jv
from scipy.fftpack import fft

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
        plt.plot(t,data[:,3])
        
        
    
def plot_gauss_data(f,fig_nr,opt='Ws',xlim=None,name=None):
    centers = np.loadtxt(f+'cv_centers.dat')
    t = centers[:,0]/1000
    c = centers[:,1]
    if np.size(centers,1)>2:
        W = centers[:,2]
    else:
        heights = np.loadtxt(f+'heights.dat')
        W = heights[:,1]
    plt.figure(fig_nr)
    plt.clf()
    if(opt=='Ws'):
        plt.plot(c,W,'x')
        plt.xlabel('$s_k$')
        plt.ylabel('$W_k~(k_\mathrm{B}T)$')
    if(opt=='st'):
        plt.plot(t,c,'x')
        plt.xlabel('$t\,(ns)$')
        plt.ylabel('$s_k$')
    if(opt=='Wt'):
        plt.plot(t,W,'x')
        plt.xlabel(r'$t\,(\mathrm{ns})$')
        plt.ylabel('$W_k$')
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
        
def plot_s_int(f, fig_nr, clear=1, opt=''):
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
        plt.plot(s,-np.log(shist))
    else:
        plt.plot(s,shist)
        plt.plot(s,Whist)
        plt.plot(s,Ehist)
    #plt.xlim([-40,40])
    ds = s[1]-s[0]
    print(sum(Whist)*ds)
    
        
    
        
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
                plt.xlabel(r'$\Gamma$')
                plt.ylabel(r'$\log\,p(\Gamma)$')            
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
    
def plot_energies(f,fig_nr,clear=1,var='beta',beta=1,linestyle='-',marker='v',label='Total energy',plot_all=False):
    data = np.loadtxt(f+'results.dat',comments='%')
    if(data.ndim==1):
        data = [data,float('NaN')*np.ones(len(data))]
    if(data.ndim>1):
        x = data[:,0]
        num_obs = round((np.size(data,1)-1)/2)
        print(num_obs)
        if(num_obs==1):
            etot = data[:,1]/en_scale
            etot_e = data[:,2]/en_scale
        if(num_obs==3):
            etot = data[:,3]/en_scale
            etot_e = data[:,4]/en_scale
        if(num_obs>3):
            etot = data[:,5]/en_scale
            etot_e = data[:,6]/en_scale
   
        plt.figure(fig_nr)
        if(clear):
            plt.clf()
        if(var=='P'):
            plt.xlabel(r'$P$')
        elif(var=='beta'):
            plt.xlabel(r'$T\,(\mathrm{K})$')
            x = 11.6045/x
        elif(var=='tau'):
            x /= beta
            #plt.grid(True)
            plt.xlabel(r'$1/\tau~(\mathrm{meV})$')
        elif(var=='tau2'):
            x = 2.0/x**2
            plt.xlabel(r'$1/\tau^2~(\mathrm{meV}^2)$')
        plt.errorbar(x,etot,etot_e,marker=marker,label=label,linestyle=linestyle)
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

        plt.ylabel(r'$\mathrm{Energy}~(\hbar\omega_0)$')
        #plt.ylabel(r'$\mathrm{Singlet~energy}~(\mathrm{meV})$')
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
        norm_shell = np.pi*((r+dr)**d-r**d)
        p /= norm_shell
        perr /= norm_shell
    return (p2,perr2)
    
    
def plot_rAB(f,fig_nr,clear=True,d=2,color='blue',marker='x',name=None,linestyle='-',show_errors=1):
    if name is None:
        name=f[-25:]
    data = np.loadtxt(f+'Pair_correlation.dat')
    r = data[:,0]
    p = data[:,1]
    p_err = data[:,2]
    (p2,perr2) = normalize(p,p_err,r,d,1)
    #p *= 1e8
    #p_err *= 1e8
    if(fig_nr>-1):
        plt.figure(fig_nr)
    if(clear):
        plt.clf()
    plt.plot(r,p2,color=color,linestyle=linestyle)
    if(show_errors):
        n = show_errors
        plt.errorbar(r[1::n],p2[1::n],perr2[1::n],linestyle='None',label=name,marker=marker,color=color)
    #if(show_errors>1):
    #    n=show_errors
    #    plt.errorbar(r[::3],p[::3],p_err[::3],linestyle='None',label=name,marker=marker,color=color)
    plt.xlabel('$r~(\mathrm{nm})$')
    plt.ylabel(r'$g(r)$')
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
            xi = -2j*x*R/a**2
            exponential = np.exp(-(2*x**2+R**2)/a**2)
            if(sym=='dis'):
            #p2part = np.exp(-0.5*(r/a)**2)
                bessels = jv(0,xi)+0.5*((2*x**2+R**2)*jv(0,xi)+1j*2*x*R*jv(1,xi))/a**2 * np.exp(-hwb)
                p2part[i]=np.absolute(sum(x*R*exponential*bessels)*dx)
                #p2part[i]=R**2*np.exp(-0.5*R**2/a**2)
                label=''
                color='r'
            if(sym=='bos'):
                p2part[i]=np.absolute(sum(x*R*exponential*jv(0,xi))*dx)
                #p2part[i]=np.absolute(sum(x*R*exponential*(2*np.pi*jv(0,2*xi)+jv(0,xi)**2))*dx)
                label=''
                color = 'b'
            if(sym=='fer'):
                p2part[i]=np.absolute(sum(x*R**3*exponential*jv(0,xi))*dx)
                label=''
                color='g'
        #color='k'
            #p2part = r**2*np.exp(-(r/a)**2)
    (p2,_)=normalize(p2part,0,r,d,1)
    dr = r[1]-r[0]
    print(p2.sum()*dr)
    #scale=1e6
    if(fig_nr>-1):
        plt.figure(fig_nr)
    plt.plot(r,p2,color=color,label=label)
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
            
def plot_1d_dist(f,fig_nr,clear=1):
    data = np.loadtxt(f+'Prob_dist1d.dat')
    r = data[:,0]
    dim = round((np.size(data,1)-1)/2)
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    for d in range(dim):
        p = data[:,2*d+1]
        perr = data[:,2*d+2]
        normalize(p,perr,r,)
        plt.errorbar(r,data[:,2*d+1],data[:,2*d+2])
    plt.xlim([-20,20])
    plt.xlabel(r'$r~(\mathrm{nm})$')
    plt.ylabel(r'$p(r)$')
    
def smooth(data, n):
    for k in range(n):
        for i in range(1,data.shape[0]-2):
            for j in range(1,data.shape[1]-2):
                data[i,j] = (data[i-1,j]+data[i,j-1]+4*data[i,j]+data[i+1,j]+data[i,j+1])*0.125
    return data

def darken(x, ):
    return x*0.8    
        
def plot_2d_dist(files,fig_nr,titles,suptitle='', stride=1,use_contour=True,
                 num_smooths=0):
    fig = plt.figure(num=fig_nr,figsize=(12,4.5))    
    #fig.set_size_inches(12,6.2)
    plt.clf()
        
    for i,f in enumerate(files):
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
        hist = smooth(histpre,num_smooths)
        norm = sum(sum(hist))*dr1*dr2
        hist = hist/norm*1000
        X,Y = np.meshgrid(r1,r2)
        Z = hist.reshape(X.shape)#.transpose()
        #print(Z.shape)
        if use_contour:

            if i==0:
                ax0 = fig.add_subplot(121)
                #Z0 = Z
            else:
                ax1 = fig.add_subplot(122,sharey=ax0)
                plt.setp(ax1.get_yticklabels(),visible=False)
                xticks = ax1.xaxis.get_major_ticks()
                #Z = (2*np.sqrt(Z)-np.sqrt(Z0))**2
                #xticks[0].label1.set_visible(False)
            #levels=np.arange(2,11,2)
            #dark_inferno=cmap_map(darken,plt.get_cmap('inferno_r'))
            plt.imshow(Z,interpolation='gaussian',extent=[r1[0],r1[-1],r2[0],r2[-1]])
            #c = plt.contour(X,Y,Z,cmap=cm.viridis)#,levels=levels)
            #plt.clabel(c, inline=1,fontsize=16,fmt='%1.1f')
            plt.xlabel(r'$x~(\mathrm{nm})$')
            if i==0:
                plt.ylabel(r'$y~(\mathrm{nm})$')
                ax0.set_title(titles[0],fontsize=28)
            else:
                ax1.set_title(titles[1],fontsize=28)
            plt.xlim([-60,60])
            plt.ylim([-30,30])
        else:
            ax = fig.add_subplot(111,projection='3d')
            surf = ax.plot_surface(X,Y,Z, rstride=stride, cstride=stride, cmap=cm.viridis)
            ax.set_zlabel(r'$p(x,y)$')
            ax.set_xlabel(r'$x~(\mathrm{nm})$')
            ax.set_ylabel(r'$y~(\mathrm{nm})$')
            ax.set_zlim([-0.05,0.5])
    #plt.set_aspect('equal')
    plt.subplots_adjust(wspace=None)
    plt.suptitle(suptitle,x=0.53,y=0.995,fontsize=32)
    if use_contour:
        ax0.set_aspect('equal')
        ax1.set_aspect('equal')
            
def free_energy_diff(f1,f2,fig_nr,beta=1.0,P=10,opt='FB',deg=1):
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
    hist_O = dE_conn[:,1]
    hist_O /= sum(hist_O)*ds
    plt.plot(s, -np.log(hist_oo), 'b', label='From oo')
    plt.plot(s, -np.log(hist_O),'k--',label='From O' )
    C = -0.041#-0.08
    #plt.plot(s, hist_oo*f_FD(s+C)*50)
    #plt.plot(s, hist_O*f_FD(-s-C)*50)
    plt.xlim([-50,50])
    plt.ylim([0,20])
    plt.xlabel(r'$s=\beta(E_O-E_{oo})$')
    plt.ylabel(r'$- \log\,p(s)$')
    #plt.legend(loc='lower right',fontsize=20)
    num=0
    den=0
    sq_oo=0
    sq_O=0
    for i,S in enumerate(s):
        num += f_FD(S+C)*hist_oo[i]
        den += f_FD(-S-C)*hist_O[i]
        sq_oo += f_FD(S+C)**2*hist_oo[i]
        sq_O += f_FD(-S-C)**2*hist_O[i]
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
    print("F_O-F_oo with f_FD:\t"+str(dF))

    #dF_p = -np.log(ds*sum(hist_oo*np.exp(-s)))
    #print("F_O-F_oo, only oo:\t"+str(dF_p))
    dF_FB = -np.log((1-np.exp(dF))/(1+np.exp(dF)))
    print("F_F-F_B with f_FD:\t"+str(dF_FB/beta))
    #dF_FBp = -np.log((1-np.exp(-dF_p))/(1+np.exp(-dF)))
    #print("F_F-F_B, only oo:\t"+str(dF_FBp))
    #print("Error in DeltaF:\t"+str(np.sqrt(err)))
    
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
        diff = b1u-b1c
        diffl = np.log(b1u)-np.log(b1c)
        plt.plot(cs,np.log(b1u),'b')
        plt.plot(cs,np.log(b1c),'r')
        #plt.plot(cs,diff,'g')
        plt.plot(cs,diffl,'c')
        poly = np.polyfit(cs,diffl,1)
        plt.plot(cs,np.polyval(poly,cs),'x--')
        C = -poly[1]/poly[0]
        plt.plot([C],[0],'bo')
        Czeros[i-1]=C
    #print(Czeros)
    #print(np.mean(Czeros))
    if(opt=='FB'):
        dF = -np.log((1-np.exp(Czeros))/(1+np.exp(Czeros))/deg)
        dF_a = -np.log((1-np.exp(np.mean(Czeros)))/(1+np.exp(np.mean(Czeros)))/deg)
    if(opt=='FD'):
        dF = -np.log((1-np.exp(Czeros))/2)
        dF_a = -np.log((1-np.exp(np.mean(Czeros)))/2)
    print("Option: "+opt)
    print("F_O-F_oo="+str(np.mean(Czeros))+"+/-"+str(np.std(Czeros)/np.sqrt(numblocks)))
    print("dF="+str(np.mean(dF/beta))+"+/-"+str(np.std(dF/beta)/np.sqrt(numblocks)))
    print("dF_alternative_mean="+str(dF_a/beta))    
    avg_sign = (1-np.exp(Czeros))/(1+np.exp(Czeros))
    print("<sign>="+str(np.mean(avg_sign))+"+/-"+str(np.std(avg_sign)/np.sqrt(numblocks)))
    
def plot_theoretical_energies(fig_nr,clear=False,d=1,hw=3.0,color='k',label=''):
    T = np.linspace(1.5,25,500)
    kB = 1/11.6045
    beta = 1.0/(kB*T)
    tau = 0.0001
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
    #Results are in units of hw
    plt.figure(fig_nr)
    if(clear):
        plt.clf()
    plt.plot(T,E,'--',color=color,label=label)
    plt.plot(T,Eb,'--',color=color)
    plt.plot(T,Ef,'--',color=color)
    
def plot_bennett_vs_T(f,fig_nr):
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
    fs=[f1,f2,f3,f4,f5]
    interac = ['noInt/','LJ/','coulomb/RW1-34/']
    beta=['beta2/']
    md=['noMetaD/','MetaD/']
    sym = ['bos/','fer/']
    P=['P10/','P20/','P40/','P50/','P60/']
    s=['s0-25/','s0-5/','s1/','s2/','s4/']
    dt=['','dt5e-3/','dt1e-3/','dt1e-4/']
    t=['t100000/']
    fa1 = f0+interac[2]+beta[0]+sym[0]+md[0]
    fa2 = f0+interac[2]+beta[0]+sym[1]+md[0]
    fa3 = f0+interac[2]+beta[0]+sym[0]+md[1]
    fa4 = f0+interac[2]+beta[0]+sym[1]+md[1]
    ff = f0+interac[2]+'tau0-04/bos/'
    flab = '../workstation_lab/metaD_test7/'
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
    f3a = '../three/distinguish/beta1/Energy_vs_P/'
    f3b = '../three/distinguish/beta2/Energy_vs_P/'

    
    fi1 ='../three/spin3half/beta1-5/woMetaD/'
    fi2 ='../coulomb170612/anisotropy3/singlet/beta1/disconnected/'
    fi3 = '../coulomb170612/anisotropy3/triplet/beta1/'
    fi4 ='../three/spin3half/beta1-5/MetaD/test/'
    fiucon = '../ideal/boson/1D/beta2/MetaD/disconnected/'
    ficonn = '../ideal/boson/1D/beta2/MetaD/connected/'
    #fiucon = '../to_send/singlet/beta1/MetaD/disconnected/'    
    #ficonn = '../to_send/singlet/beta1/MetaD/connected/'    
    #fi3 ='../ideal/distinguish/1D/beta1/'
    #fi4 ='../ideal/boson/2D/beta2/MetaD/disconnected/'
    fi5 ='../three/spin3half/beta1-5/MetaD/'
    k = 3
    fi6 ='../LJ/distinguish/'+s[k]
    fi7 ='../LJ/boson/'+s[k]
    fi8 ='../LJ/fermion/'+s[k]
    fi9 ='../ideal/boson/2D/beta1/woMetaD_longer/disconnected/'
    fi10='../ideal/fermion/2D/beta1/woMetaD_longer/'
    
    fbennett = '../ideal/fermion/1D/Energy_vs_T/FreeEnergyDiff/'

    #plot_bennett_vs_T(fbennett,1)
    #plot_theoretical_energies(1,0)
    
    #free_energy_diff(fi1,fi2,5,beta=1.34,P=20,opt='FB')

    if 0:
        plot_rAB(fi6,2,1,2,'r',marker='x',name='Disting.',linestyle='-',show_errors=1)
        plot_rAB(fi7,2,0,2,'b',marker='x',name='Boson',linestyle='-',show_errors=1)
        plot_rAB(fi8,2,0,2,'g',marker='x',name='Fermion',linestyle='-',show_errors=1)
        plt.legend(loc='upper right',fontsize=22)    
        plt.title(r'$\sigma_\mathrm{LJ}=0.5\,\mathrm{nm}$')
        plot_gauss_data(fi7,3,'Wt')
        plot_s_int(fi8,4,1,'')
        
    if 0:
        plot_2d_dist([fi2,fi3],1,[r'$\mathrm{Singlet}$',r'$\mathrm{Triplet}$'],
                     r'$\mathrm{Anisotropy}~\omega_y/\omega_x = 3$',use_contour=True,stride=10,num_smooths=0)
    
    if 1:
        #plot_cv('../sign_cv/woMetaD/',0,100000,'')
        plot_cv('../test3/',1,200000,'rew')
        #plot_energies_vs_t(f1,0,n=100000)
        plot_gauss_data(fi3,8,'Wt',name='')
        plot_gauss_data('../test2/',9,'Wt',name='test2')
        plot_s_int(fi2,2,1,'')
        #plt.ylim([-5,5])
        #plt.xlim([-1.1,1.1])
        plt.title('sing')
        plot_s_int(fi3,3,1,'')
        #plt.ylim([-5,5])
        #plt.xlim([-1.1,1.1])
        plt.title('trip')
        plot_s_int('../test/',4,1,'')
        #plt.ylim([-1,1])
        plt.xlim([-40,40])
        plt.title('test')
        plot_s_int('../test2/',5,1,'')
        #plt.ylim([-1,1])
        plt.xlim([-40,40])
        plt.title('test2')
        #Ws,pW = plot_cont(fi2 ,4,-1,100000,cv='G')
        #Ws2,pW2 = plot_cont(fi2m,5,-1,100000,cv='G')
        """plt.figure(6)
        plt.clf()
        plt.plot(Ws,pW,'b',label='Without Metadynamics',linewidth=1.5)
        plt.plot(Ws2,pW2,'r',label='With Metadynamics')
        plt.legend(loc='upper left',fontsize=22)
        plt.title('Fermion')
        plt.xlabel('$W$')
        plt.ylabel(r'$\log\,p(W)$')"""

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
