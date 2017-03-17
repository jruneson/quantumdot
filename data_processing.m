%% 01Timestep vs energy fluctuation (no thermostat)
clf; hold on;

%data01 = load('Total_energy_cl_dt0-1.dat');
e_tot1 = load('etot_dt0-1.dat');
e_tot2 = load('etot_dt0-05.dat');
e_tot3 = load('etot_dt0-02.dat');
e_tot4 = load('etot_dt0-01.dat');
e_tot5 = load('etot_dt0-005.dat');

f = @(y) log(std(y)/mean(y));

%plot(e_tot1);
%plot(e_tot2);
plot(e_tot3);
plot(e_tot4);
plot(e_tot5);
xlim([0 200])

clf;
dt = [0.1 0.05 0.02 0.01 0.005];
rel_errors = [f(e_tot1),f(e_tot2),f(e_tot3),f(e_tot4),f(e_tot5)]
plot(log(dt), rel_errors,'o-')
xlabel('log(dt)')
ylabel('log(std(E_{tot})/mean(E_{tot}))')
poly = polyfit(log(dt),rel_errors,1)
plot(dt,exp(rel_errors),'o-')

%% 02Thermostat prob distribution
clc; clf; hold on
data = load('PD_withtherm.dat');
data2 = load('PD_notherm.dat');
plot(data(:,1),data(:,2),'-');
plot(data2(:,1),data2(:,2),'-');
%xlim([-1 1])
xlabel('x')
ylabel('p(x)')
legend('With thermostat','Without thermostat','location','north')
%% 02Thermostat total energy
clf; hold on
data = load('etot_withtherm.dat');
data2 = load('etot_notherm.dat');
plot(data(:,1),data(:,2))
plot(data2(:,1),data2(:,2))
xlim([0 100])
ylabel('total energy')
xlabel('time')
title('Without thermostat. N=2, P=30, dt=0.01')

%% 03Energy estimators
clc; clf; hold on
data = load('Pot_energy.dat');
data2 = load('Kinetic_energy.dat');
data3 =load('Total_energy.dat');
plot(data(:,1),data(:,2))
plot(data2(:,1),data2(:,2))
plot(data3(:,1),data3(:,2))
plot(data(:,1),2*data(:,2),'k')
xlim([500 570])
ylim([-5 4])
legend('Pot energy','Kin energy','Tot energy', 'Tot en (with virial kin)','location','southeast')
%data(:,3) contains errors (from standard deviation)
xlabel('Time')
ylabel('Energy')
set(gca, 'FontSize',14)

%% 04Convergence in P
clf; clc; hold on
%total time 10000
P = [30 40 50 60 70 80];
dt = [0.02 0.01 0.005 0.002];
epot =     [0.490 0.511 0.515 0.507 0.523 0.601 ...
            0.526 0.512 0.569 0.508 0.511 0.556 ...
            0.569 0.510 0.547 0.523 0.563 0.533 ...
            0.520 0.549 0.548 0.518 0.540 0.558].';
epot_err = [0.062 0.104 0.078 0.109 0.108 0.151 ...
            0.096 0.106 0.159 0.081 0.110 0.195 ...
            0.128 0.070 0.192 0.179 0.249 0.224 ...
            0.123 0.194 0.153 0.162 0.202 0.198].';
ekin =     [0.537 0.506 0.533 0.527 0.477 0.592 ...
            0.558 0.486 0.601 0.501 0.491 0.518 ...
            0.500 0.476 0.376 0.467 0.529 0.479 ...
            0.513 0.542 0.507 0.618 0.510 0.470].';
ekin_err = [0.148 0.194 0.127 0.201 0.173 0.271 ...
            0.134 0.197 0.187 0.201 0.247 0.182 ...
            0.161 0.164 0.217 0.230 0.181 0.221 ...
            0.112 0.218 0.191 0.267 0.256 0.242].';
        
%total time 50000
dt = [0.02 0.01];
epot =     [0.496 0.508 0.508 0.507 0.511 0.529 ...
            0.505 0.513 0.520 0.492 0.503 0.513].';
epot_err = [0.041 0.035 0.035 0.029 0.032 0.075 ...
            0.040 0.038 0.061 0.047 0.038 0.052].';
ekin =     [0.534 0.482 0.509 0.519 0.479 0.515 ...
            0.524 0.499 0.553 0.513 0.541 0.505].';
ekin_err = [0.077 0.060 0.053 0.097 0.105 0.092 ...
            0.058 0.082 0.086 0.069 0.092 0.094].';


%% 04convergence in P and dt
clf; clc; hold on

file = importdata('run3/results.dat');
data = file.data;
P = data(:,1);
dt = data(:,2);
epot = data(:,3);
epot_err = data(:,4);
ekin = data(:,5);
ekin_err = data(:,6);
plot3d_data(P,dt,ekin,ekin_err);
zlim([0.3 0.7])
zlim([0 1])
xlabel('Number of beads')
ylabel('Timestep')
zlabel('Potential energy')
zlabel('Kinetic energy')
colormap hsv
camlight right
view(130, 30)
set(gca,'FontSize',13)

%% 05fluctuations of ekin and position
clf; hold on
data2 = load('run1/Kinetic_energy.dat');
data = load('run1/X_coord_n1p1.dat');
t = data(:,1);
x = data(:,2);
e = data2(:,2);
L = length(x);
plot(t,x,'o-')
plot(t,e,'o-')

xlim([0 15])
xlabel('Time')
ylabel('a.u.')
legend('Position of 1 bead','Kinetic energy estimator')
set(gca,'FontSize',14)
%% fft check
hold on
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = 1/(t(2)-t(1))
f = Fs*(0:(L/2))/L;
plot(f,P1);
xlim([0 3])
ylim([0 0.3])

%% 05conv in P and ttot
clf; clc; hold on
opt = 'pot';
file = importdata('run3/results.dat');
data = file.data;
P = data(:,1);
beta = data(:,2)
epot = data(:,3);
epot_err = data(:,4);
ekin = data(:,5);
ekin_err = data(:,6);
if (opt == 'pot')
    plot3d_data(P,1./beta,epot,epot_err);
    %zlim([0.4 0.6])
    zlabel('Potential energy')
elseif (opt == 'kin')
    plot3d_data(P,1./beta,ekin,ekin_err);
    zlim([-1 1])
    zlabel('Kinetic energy')
end
xlabel('Number of beads')
ylabel('T')
colormap hsv
camlight right
view(140, 10)
set(gca,'FontSize',13)
hw = 1;
E = @(kT) 0.5*hw + hw./(exp(hw./kT)-1);
T = 1./beta;
x = linspace(0,max(T));
plot3(zeros(length(x),1),x,E(x),'--')
E(max(T))

%%
clc; clf; hold on
file = importdata('run5/results.dat');
data = file.data;
%file2 = importdata('run3/results.dat');
%data2= file2.data;
beta = data(:,2);
T = 1./beta;
epot = data(:,3);
epot_err = data(:,4);
epot2 = data2(:,3);
epot_err2 = data2(:,4);
ekin = data(:,5);
ekin_err = data(:,6);
evir = data(:,7);
evir_err=data(:,8);
evir2 = data2(:,7);
evir_err2 = data2(:,8);
clf; hold on
errorbar(T,epot,epot_err,'o')
errorbar(T,ekin,ekin_err,'v')
errorbar(T,evir,evir_err,'^')
%errorbar(T,evir,evir_err,'^')
%errorbar(T,evir2,evir_err2,'o')
%plot([10 90],[0.5 0.5],'k--')
%legend('dt=T/100','dt=T/20','Target','location','southwest')
xlabel('Temperature')
ylabel('')
set(gca,'Fontsize',14)
%xlim([10 90])
%ylim([0.4 0.6])

hw = 1;
E = @(kT) 0.5*hw + hw./(exp(hw./kT)-1);
Eb = @(kT) hw*(1+1./(exp(hw./kT)-1) + 2./(exp(2*hw./kT)-1))/2;
Ef = @(kT) hw*(2+1./(exp(hw./kT)-1) + 2./(exp(2*hw./kT)-1))/2;
x = linspace(0, max(T));
plot(x,E(x),'k--')
plot(x,Eb(x),'k--')
plot(x,Ef(x),'k--')
xlim([0 3.5])
ylim([0 5])
grid on
legend('Pot','Kin','Vir', 'Theory','location','northwest')
title('Disting.part. energy')
%%
clf; clc;
data = load('run2/Pot_energy.dat');
data2= load('run2/Kin_en_virial.dat');
t = data(:,1);
pot = data(:,2);
vir = data2(:,2);
%%
clf; hold on
plot(t,pot,'o-')
plot(t,(vir-5)*100+5,'o-')
xlim([0 10])

%% 08units, 10meV
clf; hold on; clc;
type='fer';
if(type=='dis')
    folder='run4/';
elseif(type=='bos')
    folder='run5/';
elseif(type=='fer')
    folder='run6/';
else
    folder='run3/';
end

file = importdata(strcat(folder,'results.dat'));
data = file.data;
%P = [100 50 20 10 5 2].';
tau = data(:,1);%*11.6045;
beta = data(:,2);
hw=3;
data=data/hw;
epot = data(:,3);
epot_err = data(:,4);
ekin = data(:,5);
ekin_err = data(:,6);
evir = data(:,7);
evir_err = data(:,8);
x=11.6045./beta;
%x=1./tau;
errorbar(x,epot,epot_err,'o-')
errorbar(x,ekin,ekin_err,'v-')
errorbar(x,evir,evir_err,'^-')
xlim([0 max(x)])
ylim([0 6])
%plot([0 max(x)],[0.5 0.5],'k--')
%hw = 100;
k=1/11.6045;
E = @(T) 0.5*hw + hw./(exp(hw./(k*T))-1);
Eb = @(T) hw*(1+1./(exp(hw./(k*T))-1) + 2./(exp(2*hw./(k*T))-1))/2;
Ef = @(T) hw*(2+1./(exp(hw./(k*T))-1) + 2./(exp(2*hw./(k*T))-1))/2;
T = linspace(0, max(x));
if(type=='dis')
    title('Distinguishable, 2D')
elseif(type=='bos')
    title('Bosons, 2D')
elseif(type=='fer')
    title('Fermions, 2D')
end
plot(T,2*E(T)/hw,'k--')
plot(T,2*Eb(T)/hw,'k--')
plot(T,2*Ef(T)/hw,'k--')
h=legend('Potential','Kinetic','Virial','Theory','location','southeast');
set(h,'interpreter','latex');
xlabel('$T ~(\mathrm{K})$','interpreter','latex')
%xlabel('$1/\tau$ (meV)')
ylabel('Energy ($\hbar\omega_0$)','interpreter','latex')
%title('$\hbar\omega=100\,\mathrm{\mu Ha}$','interpreter','latex')
set(gca,'fontsize',18)
%title('$\hbar\omega_0=3$ meV,$\quad\beta=5$ meV$^{-1}$')
%title('$\quad\hbar\omega_0 = 3$ meV,$\quad\tau=0.15$ meV$^{-1}$')
%set(0,'defaulttextinterpreter','latex')
%set(0,'defaultaxesfontname','arial')
%ax = gca;
%xt = get(gca,'xtick');
%yt = get(gca,'YTick');
%set(ax,'XTick','fontsize',12)
xp = [0.5 0.5];
yp = [0.4 0.65];
%annotation('textarrow',xp,yp)


%%
clf; hold on; clc;
%data = load('run1/Prob_distribution.dat');
data = load('run1/X_coord_n1p1.dat');
%data = load('run1/Pot_energy_cl.dat');
%data = load('run1/Kin_energy_cl.dat');
%data = load('Pot_energy.dat');
%data = load('Kinetic_energy.dat');
%data = load('Kin_en_virial.dat');
%data = load('Pot_energy.dat')+load('Kinetic_energy.dat');
t=data(:,1);
epot = data(:,2);
plot(t,epot,'o')
xlim([0 10])
xlabel('$t$ (ps)')
ylabel('Coordinate $(a_0)$')
title('Single bead')
set(gca,'fontsize',18)
%%
hbar = 1.0545718e-34;
m = 9.10938356e-31;
a0 = 5.2917721067e-11;
Eh = 4.359744650e-18;
meV = 1.60217662e-22;
kB = 1.38064852e-23;
r = 50e-9;
%hbar*hbar/(2*m*r^2) / Eh
m/(hbar/a0/sqrt(Eh))^2;
Eh/a0^2*1e-24;
(1e-12/a0)^2;
m/(hbar/a0/sqrt(meV))^2;
meV/a0^2*1e-24;
1/(a0^2*1e24)
meV/kB
%%
folder = 'run1/';
file = importdata(strcat(folder,'results.dat'));
data = file.data;
beta = data(:,2);
epot_cl = data(:,11)/2;
epot_cl_e=data(:,12);
ekin_cl = data(:,13)/2;
ekin_cl_e=data(:,14);
clf; hold on;
T = 11.6045./beta;
errorbar(T,epot_cl,epot_cl_e,'-','linewidth',1)
errorbar(T,ekin_cl,ekin_cl_e,'-.','linewidth',1)
%plot(T,epot_cl)
%plot(T,ekin_cl)
xlabel('$T$ (K)')
ylabel('Energy (meV)')
set(gca,'fontsize',16)
title('Classical energies, single bead')
%title('Classical energies, $\tau=0.15$ meV$^{-1}$')
legend('Potential energy','Kinetic energy','location','southeast')
f = @(x) x/11.6045*0.5;
x = linspace(0,max(T));

yyaxis right
%plot(x,x/11.6045,'k--')
errorbar(T,epot_cl-f(T),epot_cl_e)
errorbar(T,ekin_cl-f(T),ekin_cl_e,'-.')
ylabel('Energy$-k_\mathrm{B}T/2$')
%ylim([-0.2 0.2])
%legend('Potential energy','Kinetic energy')