function [ h ] = plot3d_data( x,y,z,z_err )
% Plot measured data with errorbar and with interpolation surface

xlin = linspace(min(x),max(x),15);
ylin = linspace(min(y),max(y),15);
[X, Y] = meshgrid(xlin,ylin);
%u = ones(1,length(dt));
%x = kron(u,P);
%v = ones(1,length(P));
%y = kron(dt,v);
size(x)
size(y)
size(z)
f = scatteredInterpolant(x,y,z);
Z = f(X,Y);
if sum(size(Z))>0
    mesh(X,Y,Z);
end
%plot3(x,y,epot,'.','MarkerSize',15)
h = plot3_errorbar(x,y,z,z_err);
grid on; 

end

