function [h] = plot3_errorbar(x,y,z,e)
% source: http://code.izzid.com/2007/08/19/How-to-make-a-3D-plot-with-errorbars-in-matlab.html
% Jeremiah Faith 2007-08-19
% this matlab function plots 3d data using the plot3 function
% it adds vertical errorbars to each point symmetric around z
% I experimented a little with creating the standard horizontal hash
% tops the error bars in a 2d plot, but it creates a mess when you 
% rotate the plot
%
% x = xaxis, y = yaxis, z = zaxis, e = error value

% create the standard 3d scatterplot
h=plot3(x, y, z, '.k');

% looks better with large points
set(h, 'MarkerSize', 25);
hold on

% now draw the vertical errorbar for each point
for i=1:length(x)
        xV = [x(i); x(i)];
        yV = [y(i); y(i)];
        zMin = z(i) + e(i);
        zMax = z(i) - e(i);

        zV = [zMin, zMax];
        % draw vertical error bar
        h=plot3(xV, yV, zV, '-k');
        set(h, 'LineWidth', 2);
end

end

