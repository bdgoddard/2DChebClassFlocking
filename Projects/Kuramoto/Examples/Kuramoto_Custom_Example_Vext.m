% Set the external potential
optsVext.VextFn = 'VextCosSeries';
optsVext.a_n = 0.5*[1;1];  % Cosine coefficients
optsVext.phi = pi/2;
opts.optsVext = optsVext;

% Set the kernel
optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = -1;
opts.optsKernel = optsKernel;

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


