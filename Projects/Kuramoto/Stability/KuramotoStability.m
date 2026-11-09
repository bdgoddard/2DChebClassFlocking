% Housekeeping
clear all;
AddPaths;
saveDirBase = ['Kuramoto' filesep 'Stability'];

% set up the periodic interval on which to solve the eigenvalue problem
geom.N = 20;  % number of points to use for the pseudospectral method
geom.yMin = 0;    % boundaries of periodic domain
geom.yMax = 2*pi;

geomPlot = geom;   % number of points to use for the plotting
geomPlot.N = 400;

aLine = FourierLine(geom);
aLine.ComputeAll(geomPlot);

% compute the convolution integral matrix (same as in KuramotoOverdamped)
optsKernel.kernelFn = 'KuramotoKernelFourier';  % kernel function
%optsKernel.a_n = [-1;-1];  % cosine coefficients
%optsKernel.a_n = [-1;0];  % cosine coefficients
optsKernel.a_n = [-1;-1/2];  % cosine coefficients
convOpts.geom = geom;
convOpts.optsKernel = optsKernel; % note this is updated with L
convOpts.convDir = [saveDirBase filesep 'Convolutions' filesep convOpts.optsKernel.kernelFn];

Conv = ComputeConvolutionWithParameters(aLine,convOpts);

% temperature
%beta = 2.5;
beta = 0;
T = 1/beta;

% extract the relevant operators
x = aLine.Pts.y;
Dx = aLine.Diff.Dy;
DDx = aLine.Diff.DDy;
Int = aLine.Int;

% set values of r1, r2 (these are made up for now)
r1 = 0;
r2 = 0;

% get the a's
a1 = optsKernel.a_n(1);
a2 = optsKernel.a_n(2);

% exponent of the steady state and its derivative
exponent =  -a1*r1*cos(x) - a2*r2*cos(2*x);
Dexponent = a1*r1*sin(x) + 2*a2*r2*sin(2*x);

% steady state
rhoBar = exp(beta*exponent);
Z = Int*rhoBar;
rhoBar = rhoBar/Z;

% 1/rhoBar and its derivative
rhoBarInv = 1./rhoBar; % Z exp(-beta*exponent)
DrhoBarInv = -beta*Dexponent.*rhoBarInv;

% linear operator corresponding to the left hand side of the equation
L = diag(rhoBarInv)*DDx + diag(DrhoBarInv)*Dx - beta*Conv;

[evecs, evals] = eig(L);
evals = diag(evals); % output from eig is a diagonal matrix
% compute the eigenvalues
[e,sortIndex] = sort(evals);
e
% figure
% 
% % sortIndex
% % evals(sortIndex)
% 
% optsPlot.linecolor = 'r';
% optsPlot.linewidth = 2;
% 
% optsPlotExact.linecolor = 'k';
% optsPlotExact.linewidth = 2;
% optsPlotExact.linestyle = '--';
% 
% InterPol = aLine.Interp.InterPol;
% xInterp = aLine.Interp.yPlot;

% check orthogonality of eigenfunctions

% for evec = 1:size(evecs,2)
%     f = evecs(:,sortIndex(evec));
%     fInterp = InterPol*f;
%     [fMax,fMaxPos] = max(fInterp);
%     f = f/fMax;
%     shift = xInterp(fMaxPos);
%     cla
%     aLine.plot(f,optsPlot);
%     hold on
%     eval = evals(sortIndex(evec));
%     if(real(eval)<0)
%         n = round(sqrt(-eval/(2*pi)))
%         aLine.plot(cos(n*(x-shift)),optsPlotExact);
%     end
%     ylim([-1;1]);
%     pause
% end
    
