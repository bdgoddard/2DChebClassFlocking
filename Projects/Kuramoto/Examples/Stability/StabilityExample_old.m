betaPert = 0.05;

opts.geom.N = 80;
opts.geom.yMin = 0;
opts.geom.yMax = 2*pi;

data(1).optsKernel.a_n = [-1;-2];
data(1).beta_c = 1;
data(1).tMax = 100;

data(2).optsKernel.a_n = [-1;-1/2];
data(2).beta_c = 2;
data(2).tMax = 500;

data(3).optsKernel.a_n = [-1/4;-1];
data(3).beta_c = 2;
data(3).tMax = 200;

% ???
data(4).optsKernel.a_n = [-1/2;-1/4]; 
data(4).beta_c = 4;
data(4).tMax = 2000;

data(5).optsKernel.a_n = [-1;-1;-1/2];
data(5).beta_c = 2;
data(5).tMax = 200;

data(6).optsKernel.a_n = [-1;-1/4;-1];
data(6).beta_c = 2;
data(6).tMax = 200;

data(7).optsKernel.a_n = [-2;-1;-1/4];
data(7).beta_c = 1;
data(7).tMax = 200;

data(8).optsKernel.a_n = [-1/2;-1/4;-1];
data(8).beta_c = 2;
data(8).tMax = 200;


% perturbation options
optsIC.rhoICFn = 'KuramotoFourierSinePerturbationGeneral';
optsIC.kappa = 1; % hard coded
optsIC.epsilon = 0.01; % amplitude of perturbation
optsIC.shift = 0;      % phase shift of perturbation
%optsIC.shift = -pi/2; % seems to give some interesting results


nCalc = length(data);

for iCalc = 1:nCalc

%for iCalc = 4

    % kernel options
    opts.optsKernel.a_n = data(iCalc).optsKernel.a_n;

    % perturbation options
    optsIC.a = opts.optsKernel.a_n;
    optsIC.b = ones(size(opts.optsKernel.a_n));      % Fourier modes to perturb

    % set optsIC
    optsPert.optsIC = optsIC;

    % max time to run perturbed simulation
    optsPert.tMax = data(iCalc).tMax;

    % just below the critical temperature
    opts.beta = data(iCalc).beta_c - betaPert;
    optsIC.beta = opts.beta;
    KuramotoCriticalTemperature(opts,optsPert)

    % just above the critical temperature
    opts.beta = data(iCalc).beta_c + betaPert;
    optsIC.beta = opts.beta;
    KuramotoCriticalTemperature(opts,optsPert)
    
end