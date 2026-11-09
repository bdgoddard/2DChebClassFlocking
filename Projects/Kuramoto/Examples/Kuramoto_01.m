clear all

% increase the number of points to allow for an accurate simulation
geom.N = 100;  % number of points to use for the pseudospectral method
geom.yMin = 0;    % boundaries of periodic domain
geom.yMax = 1;
opts.geom = geom;

opts.geomPlot.N = 200;
opts.geomPlot.yMin = geom.yMin;
opts.geomPlot.yMax = geom.yMax;

% Set the kernel
optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = [-1];
opts.optsKernel = optsKernel;
%opts.beta = 1; % set this in the loop later instead
optsVext.VextFn = 'VextCosSeries';
optsVext.a_n = [-0.5];
% optsVext.VextFn = 'VextZero';
opts.optsVext = optsVext;

optsIC.rhoICFn = 'KuramotoSinePerturbation';  % IC function and parameters
optsIC.epsilon = 0.5;
optsIC.shift = 0;
opts.optsIC = optsIC;


opts.kappa = 1;

% Set the final time and number of timepoints to save at
%opts.toEq = false; % run for the whole time
opts.tMax = 10;
opts.nSaves = 10;

% Turn off the movie plotting
optsPlot.doPlots = true;

% do the PDE dynamics (will also show an animation of the results if you set
% optsPlot.doPlots = true)

%betaList = [0.25;0.5;1;2];
%betaList = [0.5;1;2];
%betaList = 4;
%betaList = [10];
betaList = [1;1.5;2.5;3];
nBeta = length(betaList);

for iBeta = 1:nBeta
    opts.beta = betaList(iBeta);
    dataPDE(iBeta) = KuramotoOverdamped(opts,optsPlot);
end

figure
ha = axes;

% aLine = dataPDE(1).aLine;
optsPlot.plain = true;  % turn off plotting the Chebyshev points
optsPlot.linewidth = 1.5; % set the line width
optsPlot.ha = ha; % set the axes to plot in

colorList = linspecer(nBeta);

optsPlot.linestyle = '-'; % set the line style

for iBeta = 1:nBeta
    optsPlot.linecolor = colorList(iBeta,:);  % set the line colour
    optsPlot.displayname = ['$\beta = ' num2str(betaList(iBeta)) '$'];  % set the legend name

    dataPDE(iBeta).aLine.plot(dataPDE(iBeta).rho_t(:,end),optsPlot);
    hold on
end


% fplot(@(x) exp(-cos(2*x)))

legend('interpreter','latex');
%title("W(x) = -cos(x), V(x) = - cos(2x)", 'FontWeight','Normal', "FontSize",34)