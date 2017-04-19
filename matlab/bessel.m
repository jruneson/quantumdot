%% Bessel
clc; clf; grid on
x = linspace(0,10,100);
J = besselj(0,x);
plot(x,J)
%f = 2*pi*besselj(0,2*x)+besselj(0,x).^2;
%plot(x,f)
R = 1;
J=besselj(0,-i*2*x*R);
g = R*x.*exp(-(2*x.^2+R^2)).*J;
plot(x,2*pi*g)
int = zeros(size(x));
dt = 0.05;
for t=0:dt:2*pi
    %int = int + exp(2*x*R*cos(t));
    int = int + R*x.*exp(-(2*x.^2+R^2-2*x*R*cos(t)));
end
int = int*dt;
hold on
plot(x,int)