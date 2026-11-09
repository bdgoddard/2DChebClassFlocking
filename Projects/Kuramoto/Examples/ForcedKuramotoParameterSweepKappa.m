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


% set the interaction strength to -1 (this uses the notation of the
% previous paper) - this would be 1 in the new paper, I think
opts.optsKernel.a_n = -1;

% in principal, most of the code generalises to multichromatic, but the
% exact solutions don't and it would be good to understand when/why they
% fail before moving on to the more complicated case
% opts.optsKernel.a_n = [-1;-1];

%--------------------------------------------------------------------------
% Potential and forcing options
%--------------------------------------------------------------------------

opts.h = 0;  % we're using h>0, but the code should work fine for h<0
opts.F = 0.5;

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

dKappa = 0.05;
kappaRange = (dKappa:dKappa:2.55)';

%--------------------------------------------------------------------------
% Set up the options for the Picard solver
%--------------------------------------------------------------------------

optsPicard.a = opts.optsKernel.a_n;
optsPicard.h = opts.h;
optsPicard.F = opts.F;
optsPicard.beta = opts.beta;
picardTol = 10^-9;
optsPicard.tol = picardTol;

nDim = length(opts.optsKernel.a_n);
nTests = length(kappaRange);

% store the values of rc and rs for the various kappa
kappaRC = zeros(nTests,nDim);
kappaRS = zeros(nTests,nDim);

% initial guesses - we will update these with continuation, so it really
% only matters for the first value of kappa
rc0 = 4;
rs0 = 0;

% self-consistent solve for the r's
for iTest = 1:nTests

    %kappaRange(iTest)  % uncomment this if you want to see progress (or, more likely, when it gets stuck)

    % deal with instability at the critical value if we don't have an
    % external potential
    if(~opts.h==0 || kappaRange(iTest) ~= 2) % if h = 0, skip kappa = 2

        % do some basic continuation (we know the rc and rs that solve for
        % the previous kappa, so use these as the initial guess).  In the
        % case h>0, the curves seem to be smooth so this should be a very
        % good guess.
        % DON'T DO THIS IF h=0!
        % if(iTest>1)
        %     rc0 = rc;
        %     rs0 = rs;
        % end

        % use a Picard solver to find rc and rs
        optsPicard.kappa = kappaRange(iTest);
        [rc,rs] = rForcedKuramotoPicard(rc0,rs0,aLine,optsPicard);
        kappaRC(iTest,:) = rc';
        kappaRS(iTest,:) = rs';

    end
end

disp(kappaRC(end,:))
disp(kappaRS(end,:))

%--------------------------------------------------------------------------
% Plot the results
%--------------------------------------------------------------------------


% Plot the r's against kappa

figure
plot(kappaRange,kappaRC(:,1),'-ro','DisplayName','rc1');
hold on
plot(kappaRange,kappaRS(:,1),'-mx','DisplayName','rs1');

xlabel('$\kappa$', 'interpreter','latex')
legend

% plot rho at a given kappa (set by kappaPlotPos)

figure
aRho = axes;
kappaPlotPos = length(kappaRange);
kappa = kappaRange(kappaPlotPos);
disp(['kappa = ' num2str(kappa)]);
optsRho = optsPicard;
optsRho.kappa = kappa;
rhoPlot = ForcedKuramotoRhoFromR(kappaRC(kappaPlotPos,1),kappaRS(kappaPlotPos,1),aLine,optsRho);

optsPlot.plain = true;
optsPlot.linewidth = 2;
optsPlot.linestyle = '-';
optsPlot.color = 'r';
aLine.plot(rhoPlot,optsPlot);


% do the parameter passing trick for fsolve to find a, as in Greg's paper
beta = opts.beta;
eta = opts.h;
% we set kappa above already
f = @(x)aRHS(x,beta,eta,kappa);

% plot the two parts of the self-consistent equation to visualise the roots
x = -4:0.01:4;
figure
plot(x,f(x),'r')
hold on
plot(x,eta + kappa * besseli(1,x)./besseli(0,x),'b')

% note that we need to take the initial guess for a to be suitably large
% so that we don't just find the a=0 solution
% We take the sign to be the same sign as h.

optsFSolve = optimoptions(@fsolve,'FunctionTolerance',1e-10);

% deal with the fact that sign(0) = 0
if(eta==0)
    a0 = 4;
else
    a0 = 4*sign(eta);
end
aMin = fsolve(f,a0,optsFSolve);
f(aMin) % this should be small

% construct the exact equilibrium solution and plot it

q = aLine.Pts.y;
Int = aLine.Int;
rhoInf = exp(aMin*cos(q));
Z = Int*rhoInf;
rhoInf = rhoInf/Z;

axes(aRho);
hold on
optsPlot.plain = true;
optsPlot.linewidth = 2;
optsPlot.linestyle = '--';
optsPlot.color = 'b';
aLine.plot(rhoInf,optsPlot);


function F = aRHS(a,beta,eta,kappa)
    % https://doi.org/10.1007/s00205-021-01648-1
    % Lemma 1.26
    % self-consistent equation for a is given by (4.2), where we introduce
    % \kappa (the paper only uses \beta)
    F = a - beta * (eta + kappa * besseli(1,a)/besseli(0,a));
end
