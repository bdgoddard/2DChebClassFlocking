% Set the external potential
eta = 0.5;
optsVext.VextFn = 'VextCosSeries';
optsVext.a_n = [0;-eta];  % Cosine coefficients
optsVext.phi = 0;
opts.optsVext = optsVext;

% Set the kernel
optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = -1;
opts.optsKernel = optsKernel;

opts.beta = 0.1;

% Set the final time and number of timepoints to save at
opts.tMax = 5;
opts.nSaves = 40;

% do the PDE dynamics (will also show an animation of the results)
dataPDE = KuramotoOverdamped(opts);

% plot the external potential
aLine = dataPDE.aLine;
theta = aLine.Pts.y;
VextFn = str2func(opts.optsVext.VextFn);
V = VextFn(theta,opts.optsVext);
optsPlot.plain = true;
figure
aLine.plot(V,optsPlot);


