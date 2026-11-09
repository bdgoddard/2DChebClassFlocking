% Set the kernel
optsKernel.kernelFn = 'KuramotoKernelCosSquared';
optsKernel.a_n = -1;
opts.optsKernel = optsKernel;

% Set the final time and number of timepoints to save at
opts.tMax = 30;
opts.nSaves = 40;

% Turn off the movie plotting
optsPlot.doPlots = false;

% do the PDE dynamics (will also show an animation of the results if you set
% optsPlot.doPlots = true)
dataPDECosSq = KuramotoOverdamped(opts,optsPlot);

% this should be equivalent as the constant term in 
% cos^2(q) = 1/2(1 + cos(2*q) doesn't matter in the dynamics
optsKernel.kernelFn = 'KuramotoKernelFourier';  % kernel function
optsKernel.a_n = [0;-1/2];                      % cosine coefficients
opts.optsKernel = optsKernel;

dataPDECos2 = KuramotoOverdamped(opts,optsPlot);

rho1 = dataPDECosSq.rho_t;
rho2 = dataPDECos2.rho_t;
err = max(abs(rho1(:) - rho2(:)));

disp(['L^\infty error between the two dynamics: ' num2str(err)]);

figure
ha = axes;

aLine = dataPDECosSq.aLine;
optsPlot.plain = true;  % turn off plotting the Chebyshev points
optsPlot.linewidth = 1.5; % set the line width
optsPlot.ha = ha; % set the axes to plot in

colorList = linspecer(2);

optsPlot.linecolor = colorList(1,:);  % set the line colour
optsPlot.displayname = '$\cos^2(x)$';  % set the legend name
optsPlot.linestyle = '-'; % set the line style

aLine.plot(rho1(:,end),optsPlot);

hold on

optsPlot.linecolor = colorList(2,:);  % set the line colour
optsPlot.displayname = '$\cos(2x)$';  % set the legend name
optsPlot.linestyle = '--'; % set the line style

aLine.plot(rho2(:,end),optsPlot);
legend('interpreter','latex');