betaPert = 0.02;

opts.geom.N = 80;
opts.geom.yMin = 0;
opts.geom.yMax = 2*pi;

data(1).optsKernel.a_n = [-1;-2];
data(1).beta_c = 1;
data(1).tMax = 100;

data(2).optsKernel.a_n = [-1;-1/2];
data(2).beta_c = 2;
data(2).tMax = 1500;

data(3).optsKernel.a_n = [-1/4;-1];
data(3).beta_c = 2;
data(3).tMax = 400;

data(4).optsKernel.a_n = [-1/2;-1/4]; 
data(4).beta_c = 4;
data(4).tMax = 4000;

data(5).optsKernel.a_n = [-1;-1;-1/2];
data(5).beta_c = 2;
data(5).tMax = 500;

data(6).optsKernel.a_n = [-1;-1/4;-1];
data(6).beta_c = 2;
data(6).tMax = 1000;

data(7).optsKernel.a_n = [-2;-1;-1/4];
data(7).beta_c = 1;
data(7).tMax = 500;

data(8).optsKernel.a_n = [-1/2;0;0;0;0;-1/6];
data(8).beta_c = 4;
data(8).tMax = 5000;

data(9).optsKernel.a_n = [-1;-1/2;-1/3;-1/4];
data(9).beta_c = 2;
data(9).tMax = 2000;

data(10).optsKernel.a_n = [-1/4;-1/2;-1/6;-1];
data(10).beta_c = 2;
data(10).tMax = 200;

data(11).optsKernel.a_n = [-1/4;-1/10;-1/5;-1/6;-1/2];
data(11).beta_c = 4;
data(11).tMax = 200;

data(12).optsKernel.a_n = [-3;-2];
data(12).beta_c = 2/3;
data(12).tMax = 200;

%---------------------
% There are the ones in Benedetta's table:
%---------------------

data(13).optsKernel.a_n = [0;0;0;-1/4;0;-1/6];
data(13).beta_c = 8;
data(13).tMax = 1000;

data(14).optsKernel.a_n = [0;0;0;-3;0;-2];
data(14).beta_c = 2/3;
data(14).tMax = 20;

data(15).optsKernel.a_n = [0;0;-1/8;-1/2];
data(15).beta_c = 4;
data(15).tMax = 400;

data(16).optsKernel.a_n = [-1/4;0;-1/6;0;-1/2];
data(16).beta_c = 4;
data(16).tMax = 200;

data(17).optsKernel.a_n = [0;0;0;-1;0;-1/2];
data(17).beta_c = 2;
data(17).tMax = 200;

data(18).optsKernel.a_n = [-1/10;0;0;-1/4;0;-1/6];
data(18).beta_c = 8;
data(18).tMax = 2000;

data(19).optsKernel.a_n = [-1/5;0;0;-1/4;0;-1/6];
data(19).beta_c = 8;
data(19).tMax = 1500;

%---------------------

% perturbation options
optsIC.rhoICFn = 'KuramotoFourierSinePerturbationGeneral';
optsIC.kappa = 1; % hard coded
optsIC.epsilon = 0.01; % amplitude of perturbation
optsIC.shift = 0;      % phase shift of perturbation
%optsIC.shift = -pi/2; % seems to give some interesting results

optsPlot = struct;
optsPlot.doPlots = true;

nCalc = length(data);

for iCalc = 1:nCalc
    
    % kernel options
    opts.optsKernel.a_n = data(iCalc).optsKernel.a_n;

    % perturbation options
    optsIC.a = opts.optsKernel.a_n;
    optsIC.b = ones(size(opts.optsKernel.a_n));      % Fourier modes to perturb
    optsIC.r = zeros(size(optsIC.a)); % test the constant solution

    % set optsIC
    opts.optsIC = optsIC;

    % max time to run perturbed simulation
    opts.tMax = data(iCalc).tMax;

    % just below the critical temperature (should go back to constant)
    opts.beta = data(iCalc).beta_c - betaPert;
    opts.optsIC.beta = opts.beta;
    dataOut(iCalc,1) = KuramotoOverdamped(opts,optsPlot);

    % just above the critical temperature (should go to peaked state)
    opts.beta = data(iCalc).beta_c + betaPert;
    opts.optsIC.beta = opts.beta;
    dataOut(iCalc,2) = KuramotoOverdamped(opts,optsPlot);
    
end