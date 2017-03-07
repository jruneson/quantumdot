%% Timestep vs energy fluctuation (no thermostat)
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

%% Thermostat prob distribution
clc; clf; hold on
data = load('PD_withtherm.dat');
data2 = load('PD_notherm.dat');
plot(data(:,1),data(:,2),'-');
plot(data2(:,1),data2(:,2),'-');
%xlim([-1 1])
xlabel('x')
ylabel('p(x)')
legend('With thermostat','Without thermostat','location','north')
%% Thermostat total energy
clf; hold on
data = load('etot_withtherm.dat');
data2 = load('etot_notherm.dat');
plot(data(:,1),data(:,2))
plot(data2(:,1),data2(:,2))
xlim([0 100])
ylabel('total energy')
xlabel('time')
title('Without thermostat. N=2, P=30, dt=0.01')

%% 
clc; clf; hold on
%data = load('Prob_distribution.dat');
data = load('X_coordinate_n1p1.dat');
plot(data(:,1),data(:,2))
%data(:,3) contains errors (from standard deviation)
