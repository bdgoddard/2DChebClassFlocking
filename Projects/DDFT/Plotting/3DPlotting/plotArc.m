function P = plotArc(a,b,x0,y0,r)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians, 
% b is end of arc in radians, 
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plotArc(pi/4,3*pi/4,9,-4,3)
t = linspace(a,b,10);
x = r*cos(t) + x0;
y = r*sin(t) + y0;
x = [x x0 x(1)];
y = [y y0 y(1)];
P = fill(x,y,'r');
%axis([x0-r-1 x0+r+1 y0-r-1 y0+r+1]) 
%axis equal;
end