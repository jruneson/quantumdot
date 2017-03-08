function [ h ] = plot3d_data( P,dt,z,z_err )
% Plot measured data with errorbar and with interpolation surface

xlin = linspace(min(P),max(P),15);
ylin = linspace(min(dt),max(dt),15);
[X, Y] = meshgrid(xlin,ylin);
%u = ones(1,length(dt));
%x = kron(u,P);
%v = ones(1,length(P));
%y = kron(dt,v);
x = P;
y = dt;
size(x)
size(y)
size(z)
f = scatteredInterpolant(x,y,z);
Z = f(X,Y);
mesh(X,Y,Z)
%plot3(x,y,epot,'.','MarkerSize',15)
h = plot3_errorbar(x,y,z,z_err)
grid on; 

end

