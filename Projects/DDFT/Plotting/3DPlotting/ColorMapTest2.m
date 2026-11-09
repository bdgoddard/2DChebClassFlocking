close all
clear all


% Generate some surface data.
[X1,Y1,Z1] = peaks(30);
[X2,Y2,Z2] = peaks(20);
Z2 = Z2 + 4*rand(size(Z2));
% Produce the two surface plots.
h(1) = surf(X1,Y1,Z1);
hold on
h(2) = pcolor(X2,Y2,Z2);
hold off

% Move the pcolor to Z = -10.
% The 0*Z is in the statement below to insure that the size
% of the ZData does not change.
set(h(2),'ZData',-10 + 0*Z2)
set(h(2),'FaceColor','interp','EdgeColor','interp')
view(3)


% Scale the CData (Color Data) of each plot so that the 
% plots have contiguous, nonoverlapping values.  The range 
% of each CData should be equal. Here the CDatas are mapped 
% to integer values so that they are easier to manage; 
% however, this is not necessary.
% Initially, both CDatas are equal to Z.
m1 = 64;  % number of elements in each color map
m2 = 10;   

% Define a colormap that uses the cool colormap and 
% the gray colormap and assign it as the Figure's colormap.
colormap([cool(m1);gray(m2)])


cmin1 = min(Z1(:));
cmax1 = max(Z1(:));
cmin2 = min(Z2(:));
cmax2 = max(Z2(:));

% CData for surface
C1 = min(m1,round((m1-1)*(Z1-cmin1)/(cmax1-cmin1))+1); 
C2 = m1 + min(m2,round((m2-1)*(Z2-cmin2)/(cmax2-cmin2))+1); 
% CData for pcolor
%C2 = 64+C1;
% Update the CDatas for each object.
set(h(1),'CData',C1);
set(h(2),'CData',C2);
% Change the CLim property of axes so that it spans the 
% CDatas of both objects.
caxis([min(C1(:)) max(C2(:))])