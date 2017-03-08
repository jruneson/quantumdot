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


%%
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