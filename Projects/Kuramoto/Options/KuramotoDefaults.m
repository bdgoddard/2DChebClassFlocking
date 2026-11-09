function default = KuramotoDefaults()

%--------------------------------------------------------------------------
% Domain and PDE domain setup
%--------------------------------------------------------------------------

geom.N = 50;  % number of points to use for the pseudospectral method
geom.yMin = 0;    % boundaries of periodic domain
geom.yMax = 2*pi;

geomPlot = geom;   % number of points to use for the plotting
geomPlot.N = 400;

default.geom = geom;
default.geomPlot = geomPlot;

%--------------------------------------------------------------------------
% Whether to use the log form of the PDE solver
%--------------------------------------------------------------------------

default.doLog = true;

%--------------------------------------------------------------------------
% Final time and number of save points
%--------------------------------------------------------------------------

default.tMax = 200;
default.nSaves = 20;

%--------------------------------------------------------------------------
% Tolerances for the PDE solver
%--------------------------------------------------------------------------

default.relTol = 10^-9;  % relative and absolute tolerances
default.absTol = 10^-9;

default.toEq = true;     % determine if stop once the rhs is small
default.eqTol = 10^-7;   % tolerance to stop as we're 'at equilibrium'

%--------------------------------------------------------------------------
% Physical parameters
%--------------------------------------------------------------------------

default.beta = 1;  % temperature
default.kappa = 1; % interaction strength

%--------------------------------------------------------------------------
% Kernel parameters
%--------------------------------------------------------------------------

optsKernel.kernelFn = 'KuramotoKernelFourier';  % kernel function
optsKernel.a_n = [-1];                       % cosine coefficients
default.optsKernel = optsKernel;

%--------------------------------------------------------------------------
% Initial condition parameters
%--------------------------------------------------------------------------

optsIC.rhoICFn = 'KuramotoSinePerturbation';  % IC function and parameters
optsIC.epsilon = 0.01;
optsIC.shift = -pi/2;
default.optsIC = optsIC;

%--------------------------------------------------------------------------
% External potential parameters
%--------------------------------------------------------------------------

% optsVext.VextFn = 'VextZero';
% default.optsVext = optsVext;

optsVext.VextFn = 'VextCos';
optsVext.alpha = 0;
optsVext.phi = pi/2;
default.optsVext = optsVext;

%--------------------------------------------------------------------------
% Stochastic numerics parameters
%--------------------------------------------------------------------------

%optsSDE.nParticles = 500;
optsSDE.nParticles = 200;
optsSDE.nSteps = default.tMax*20;
optsSDE.nRuns = 1;
optsSDE.nSamples = optsSDE.nRuns;
optsSDE.thin = 1;
optsSDE.burnin = 100;
default.optsSDE = optsSDE;

end