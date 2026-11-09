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
% opts.optsKernel.a_n = [-1;-1];
opts.optsKernel.a_n = [-1;-1/2];
%  opts.optsKernel.a_n = [-1;0];
%opts.optsKernel.a_n = [1;-1];
% opts.optsKernel.a_n = [0;-1];


%--------------------------------------------------------------------------
% Find how the r's change (starting from an initial guess) as we change
% kappa
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Set up the domain
%--------------------------------------------------------------------------

aLine = FourierLine(opts.geom);
aLine.ComputeAll(opts.geomPlot);

%--------------------------------------------------------------------------
% Set up the range of kappa
%--------------------------------------------------------------------------

dKappa = 0.01;
kappaRange = (dKappa:dKappa:4)';

%--------------------------------------------------------------------------
% Set up the options for the Picard solver
%--------------------------------------------------------------------------

optsPicard.a = opts.optsKernel.a_n;
optsPicard.beta = opts.beta;
picardTol = 10^-9;
optsPicard.tol = picardTol;

nDim = length(opts.optsKernel.a_n);
nTests = length(kappaRange);
kappaTests = zeros(nTests,nDim);

% initial guess -- this will affect which set of r's we converge to
r0 = [0.75;0.25];

% self-consistent solve for the r's
for iTest = 1:nTests
    optsPicard.kappa = kappaRange(iTest);
    kappaTests(iTest,:) = rKuramotoPicard(r0,aLine,optsPicard)';
end
    
% The method is a bit unstable at the critical kappa -- this will depend
% very much on the choices of the a's
% if(any(kappaRange==2))
%     kappaTests(kappaRange==2,:) = [0;0];
% end

%--------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------

plot(kappaRange,kappaTests(:,1),'-o','DisplayName','r1');
hold on
plot(kappaRange,kappaTests(:,2),'-o','DisplayName','r2');
legend
