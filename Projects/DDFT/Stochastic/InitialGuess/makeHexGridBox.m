function r = makeHexGridBox(optsPhys)
% put particles on a hexagonal grid
% used to get initial configuration for sampling

nParticles = optsPhys.nParticlesS; % only works for a single species so this should be a scalar

sigma = optsPhys.sigma;     % particle diameter
R = sigma/2;

xMin = 0 + R;           % minimum centre position in x
xMax = optsPhys.L1 - R; % maximum centre position in x
yMin = 0 + R;           % minimum centre position in y
yMax = optsPhys.L2 - R; % maximum centre position in y

x1 = (xMin : sigma : xMax);     % points separated by sigma for odd rows
x2 = (xMin + R : sigma : xMax); % and shifted by R for even rows

y1 = yMin : sqrt(3)*sigma : yMax; % points separated by sqrt(3)*sigma for odd rows
y2 = yMin + sqrt(3)*R : sqrt(3)*sigma : yMax; % and shifted by sqrt(3)*R for even rows

% generate points with meshgrid for odd and even rows and then combine
[X1,Y1] = meshgrid(x1,y1);
[X2,Y2] = meshgrid(x2,y2);

% sort the list according to r1
r1 = [X1(:);X2(:)];
r2 = [Y1(:);Y2(:)];
[r1,sortOrder] = sort(r1);
r2 = r2(sortOrder);

% create the full list
r = [r1,r2];


nGridpoints = length(r(:,1));
if(nGridpoints < nParticles)
    error('makeHexGridBox: too many particles to fit in box');
else
    skipRate = floor(nGridpoints/nParticles);
    toRemoveList = 1 : skipRate : nGridpoints;
    r = r(toRemoveList,:);
    r = r(1:nParticles,:);
end


end

