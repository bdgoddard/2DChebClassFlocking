% You can overwrite any of the options in 
% -- KuramotoDefaults
% or
% -- KuramotoDefaultsPhaseSpace
% by setting something like
% optsKernel.kernelFn = 'KuramotoKernelFourier';
% optsKernel.a_n = [-1];
% opts.kappa = 5;
%
% Note that if you want to overwrite something in a structure, such as
% optsKernel, optsIC, or optsVext, you need to specify the whole structure

clear all
AddPaths;

% plotting options
optsPlot.doMp4 = false;     % make an mp4 movie
optsPlot.doSnapshots = true;  % plot snapshots
optsPlot.doPlots = false;  % show/hide plotting within the PDE solver

% Do the overdamped dynamics
opts.phasespace = false;
optsPlot.mp4File = 'example_overdamped.mp4';

opts.tMax = 20;
opts.nSaves = 8;

optsSDE.nParticles = 500;
optsSDE.nSteps = opts.tMax*100;
optsSDE.nRuns = 1;
optsSDE.nSamples = optsSDE.nRuns;
optsSDE.thin = 1;
optsSDE.burnin = 100;
opts.optsSDE = optsSDE;

% Set the external potential
optsVext.VextFn = 'VextCosSeries';
optsVext.a_n = 0.5*[1;1];  % Cosine coefficients
optsVext.phi = pi/2;
opts.optsVext = optsVext;

% Set the kernel
optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = -1;
opts.optsKernel = optsKernel;

[dataPDE_over,dataSDE_over] = KuramotoComparison(opts,optsPlot);

% do the underdamped dynamics
% opts.phasespace = true;
% optsPlot.mp4File = 'example_underdamped.mp4';
% 
% [dataPDE_under,dataSDE_under] = KuramotoComparison(opts,optsPlot);