clear all
close all
AddPaths;

% plotting options
optsPlot.doMp4 = false;     % make an mp4 movie
optsPlot.doSnapshots = true;  % plot snapshots
optsPlot.doPlots = false;  % show/hide plotting within the PDE solver

% Do the overdamped dynamics
opts.phasespace = false;

% time options
opts.tMax = 2;
opts.nSaves = 4;

% SDE options
optsSDE.nParticles = 2000;
optsSDE.nSteps = opts.tMax*100;
optsSDE.nRuns = 1;
optsSDE.nSamples = optsSDE.nRuns;
optsSDE.thin = 1;
optsSDE.burnin = 100;
opts.optsSDE = optsSDE;

% Set the external potential
optsVext.VextFn = 'VextCosSeries';
optsVext.a_n = 0.5;  % Cosine coefficients
optsVext.phi = pi/2;
opts.optsVext = optsVext;

% Set the kernel
optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = -1;
opts.optsKernel = optsKernel;

% number of bins for histogramming the SDE
optsPlot.nBins = 10;

% file to save the snapshots in
optsPlot.doPDF = true;
optsPlot.snapshotsFile = 'snapshotExample.pdf';

[dataPDE,dataSDE] = KuramotoComparison(opts,optsPlot);
