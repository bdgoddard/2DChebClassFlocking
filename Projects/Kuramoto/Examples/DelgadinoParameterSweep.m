%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------

clear all
% Ensure we have dirData as a variable
AddPaths;
% Set up save location
saveDirBase = [filesep 'Kuramoto' filesep 'Overdamped'];

%--------------------------------------------------------------------------
% Load some sensible defaults
%--------------------------------------------------------------------------

opts = KuramotoDefaults();

%--------------------------------------------------------------------------
% Kernel options
%--------------------------------------------------------------------------
%
% Note that this should work for a general kernel, of arbitrary length.  We
% are setting the fourier modes of \sum a_n cos(nq) in a vector
% [a_1; a_2; ...; a_n]

% some examples
%  opts.optsKernel.a_n = [-1;0];
opts.optsKernel.a_n = [-1;-1];
%  opts.optsKernel.a_n = [-1;1];
%opts.optsKernel.a_n = [1;-1];
% opts.optsKernel.a_n = [-1;1;-1];

% interaction strength
%opts.kappa = 2.05;
%opts.kappa = 3;
opts.kappa = 4;
%opts.kappa = 5;
%opts.kappa = 7;

%--------------------------------------------------------------------------
% Find the possible sets of r that solve the self-consistency equations
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Set up the domain
%--------------------------------------------------------------------------

aLine = FourierLine(opts.geom);
aLine.ComputeAll(opts.geomPlot);

q = aLine.Pts.y;
Int = aLine.Int;

%--------------------------------------------------------------------------
% Compute the convolutions for a Picard solver (same as used in for the PDE)
%--------------------------------------------------------------------------

% need this for the free energy (not for the Picard solvers as the
% convolutions have already been computed explicitly)
convOpts.geom = opts.geom;
convOpts.optsKernel = opts.optsKernel;
convOpts.convDir = [saveDirBase filesep 'Convolutions' filesep convOpts.optsKernel.kernelFn];
aLine.Conv = ComputeConvolutionWithParameters(aLine,convOpts);

%--------------------------------------------------------------------------
% Set up a grid of initial guesses for r, which we hope will lead to all
% possible solutions of the self-consistency equations
%--------------------------------------------------------------------------

% the 1D grid for each r (we could generalize this to different values, but
% this seems good enough)
rGrid = -1:0.2:1;

% replicated the grid for all dimensions
% warning: obscure matlab code follows!
nDim = length(opts.optsKernel.a_n);
ndgridInput = cell(1,nDim);
for iDim = 1:nDim
    ndgridInput{iDim} = rGrid;
end
ndgridOutput = cell(size(ndgridInput));
[ndgridOutput{:}] = ndgrid(ndgridInput{:});

rGridFull = zeros(length(ndgridOutput{1}(:)),nDim);

for iDim = 1:nDim
    rGridFull(:,iDim) = ndgridOutput{iDim}(:);
end

%--------------------------------------------------------------------------
% Use the picard solver starting from each element of rGridFull
%--------------------------------------------------------------------------

picardTol = 10^-9;

% copy some options over to the options sent to the Picard solver
optsPicard.a = opts.optsKernel.a_n;
optsPicard.beta = opts.beta;
optsPicard.kappa = opts.kappa;
optsPicard.tol = picardTol;

% run the Picard solver
nTests = size(rGridFull,1);
rTests = zeros(size(rGridFull));
for iTest = 1:nTests
    rTests(iTest,:) = rKuramotoPicard(rGridFull(iTest,:)',aLine,optsPicard)';
end

%--------------------------------------------------------------------------
% Clean up the results
%--------------------------------------------------------------------------

% make r1<=0; the signs of all odd numbered coefficients can all switch,
% which just results in a function shifted by pi
r1PosMask = ( rTests(:,1) > 0 );
rTests(r1PosMask,1:2:nDim) = - rTests(r1PosMask,1:2:nDim);

% sort out zero entries; small numbers mess up uniquetol as it's done by a
% relative error
rTests(abs(rTests) < 100*picardTol) = 0;

% find unique entries up to a tolerance
uTol = 10^-5;
[rUnique,rIndex] = uniquetol(rTests,uTol,'ByRows',true);

% display the results
rUnique

% plot the set of unique rho's
nUnique = size(rUnique,1);
a_n = opts.optsKernel.a_n;

figure('Position',[100,100,1200,800]);
optsPlot.plain = true;
optsPlot.linewidth = 1.5;

optsEq.a = opts.optsKernel.a_n;
optsEq.beta = opts.beta;
optsEq.kappa = opts.kappa;

for iUnique = 1:nUnique
    r_i = rUnique(iUnique,:)';
    optsEq.r = r_i;
    displayname = [];
    for iDim = 1:nDim
        displayname = cat(2,displayname,['r' num2str(iDim) ' = ' num2str(r_i(iDim))]);
        if(iDim<nDim)
            displayname = cat(2,displayname,', ');
        end
    end
    optsPlot.displayname = displayname;
    rho = KuramotoRhoFourier(q,optsEq);
    rho = rho/(Int*rho);
    
    % compute the steady state starting from the constructed rho via a
    % Picard solver (this should be the same as that constructed above, but
    % this is a good sanity check)
    aLine.plot(rho,optsPlot);
    test_r_i = rKuramotoPicard(r_i,aLine,optsPicard)';
    hold on    
end    
legend


%--------------------------------------------------------------------------
% Perturb the various fixed point rho's and run the PDE to check for
% stability
%--------------------------------------------------------------------------

% clear the initial condition from the options (set by default)
opts = rmfield(opts,'optsIC');

% construct an initial condition, which takes the exact rho, with the
% computed r values, and perturbs it (multiplicatively) with a sine wave
optsIC.rhoICFn = 'KuramotoFourierSinePerturbation';
optsIC.a = opts.optsKernel.a_n;
optsIC.kappa = opts.kappa;
optsIC.beta = opts.beta;
optsIC.epsilon = 0.01; % amplitude of perturbation
optsIC.shift = 0;      % phase shift of perturbation
%optsIC.shift = -pi/2; % seems to give some interesting results
opts.optsIC = optsIC; % append initial conditions to options

% set the final time - this is highly dependent on the choice of parameters
opts.tMax = 250; % kappa = 4, a_n = [-1,-1]
%opts.tMax = 5000; % kappa = 8, a_n = [-1;-1]
%opts.tMax = 500; % kappa = 4, a_n = [-1,-1]
%opts.tMax = 10; % kappa = 5, a_n = [-1,1]
% opts.tMax = 10; % kappa = 7, a_n = [-1,1]

% set the number of save points (to make a nice movie and plots of the free energy)
opts.nSaves = 50;

% turn off the plotting
optsPlot.doPlots = false;

% run the PDE dynamics using each value of r
for iUnique = 1:nUnique
    opts.optsIC.r = rUnique(iUnique,:)';
    data(iUnique) = KuramotoOverdamped(opts,optsPlot);
end

% plot the free energies in the same colours as the rho plots above
hf = figure('Position',[100,100,1200,800]);
maxF = -inf;
minF = inf;

for iUnique = 1:nUnique
    figure(hf);
    displayname = [];
    for iDim = 1:nDim
        displayname = cat(2,displayname,['r' num2str(iDim) ' = ' num2str(rUnique(iUnique,iDim))]);
        if(iDim<nDim)
            displayname = cat(2,displayname,', ');
        end
    end

    plot(data(iUnique).plotTimes,data(iUnique).F,'linewidth',1.5,'DisplayName',displayname);
    maxF = max(maxF,max(data(iUnique).F));
    minF = min(minF,min(data(iUnique).F));
    hold on
end

ylim([minF,maxF]);
legend

% seems that the two single clusters are stable, the other two are unstable
% look at sin(q-pi/2) perturbation - this seems to keep the two-cluster
% state stable (presumably because there's some symmetry in ther
% perturbation)


