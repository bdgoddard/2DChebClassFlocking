%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------

% clear all
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

%opts.optsKernel.a_n = [-1;-1/2];
opts.optsKernel.a_n = [-1;-1];
opts.beta = 1;

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

%kappaList = [1.5:0.01:2.5];
kappaList = [1.5:0.1:2.5];
%kappaList = 2.5;
nKappa = length(kappaList);

data = struct;

optsEq.a = opts.optsKernel.a_n;
optsEq.beta = opts.beta;


for iKappa = 1:nKappa
    opts.kappa = kappaList(iKappa);

    if(opts.kappa == 2)
        data(iKappa).kappa = opts.kappa;
        data(iKappa).rValues = [];
        data(iKappa).rhoList = [];
        continue;
    end
    
    %--------------------------------------------------------------------------
    % Find the possible sets of r that solve the self-consistency equations
    %--------------------------------------------------------------------------
    
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

    picardTol = 10^-10;

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
    rTests(abs(rTests) < 1000*picardTol) = 0;

    % find unique entries up to a tolerance
    uTol = 10^-5;
    [rUnique,~] = uniquetol(rTests,uTol,'ByRows',true);
    
    rUnique = sortrows(rUnique,'descend'); % sort by r1
    
    
    % display the results
    rUnique
    data(iKappa).kappa = opts.kappa;
    data(iKappa).rValues = rUnique;    
    
    nUnique = size(rUnique,1);
    rhoList = zeros(length(q),nUnique);
   
    for iRho = 1:nUnique
        optsEq.r = rUnique(iRho,:)';
        optsEq.kappa = opts.kappa;
        rho = KuramotoRhoFourier(q,optsEq);
        rho = rho/(Int*rho);
        rhoList(:,iRho) = rho;
    end
    data(iKappa).rhoList = rhoList;
    
end

% clear the initial condition from the options (set by default)
opts = rmfield(opts,'optsIC');

% construct an initial condition, which takes the exact rho, with the
% computed r values, and perturbs it (multiplicatively) with a sine wave
%optsIC.rhoICFn = 'KuramotoFourierSinePerturbation';
optsIC.rhoICFn = 'KuramotoFourierSinePerturbationGeneral';
 
optsIC.a = opts.optsKernel.a_n;
optsIC.beta = opts.beta;
optsIC.b = [1;1];      % strength of perturbation in different frequencies
optsIC.epsilon = 0.01; % amplitude of perturbation
optsIC.shift = 0;      % phase shift of perturbation
%optsIC.shift = -pi/2; % seems to give some interesting results
opts.optsIC = optsIC;  % append initial conditions to options

% set the final time - this is highly dependent on the choice of parameters
opts.tMax = 1000;
% set the number of save points (to make a nice movie and plots of the free energy)
opts.nSaves = 50;

% turn off the plotting
optsPlot.doPlots = false;

if(optsPlot.doPlots)
    figure
end

% run the PDE dynamics using each value of r

figure('Position',[100,100,1200,800],'DefaultAxesFontSize',20);
colourList = {'r','b','m'};
markerList = {'o','x','*'};

stableTol = 10^-2;

for iKappa = 1:nKappa
    
    opts.kappa = data(iKappa).kappa;
    opts.optsIC.kappa = opts.kappa;
    rValues = data(iKappa).rValues;
    for iUnique = 1:size(rValues,1)

        opts.optsIC.r = rValues(iUnique,:)';
        
%         opts.optsIC.r
        
        dataPDETemp = KuramotoOverdamped(opts,optsPlot);
        pause(1) % ensure different file names
        
        dataPDETemp.rhoPicard = data(iKappa).rhoList(:,iUnique);
        dataPDE(iKappa,iUnique) = dataPDETemp;

        rhoDiff = max(abs(dataPDE(iKappa,iUnique).rhoPicard - dataPDE(iKappa,iUnique).rho_t(:,end)));
        
        stable(iKappa,iUnique) = rhoDiff < stableTol;
        
        for iR = 2:2
        %for iR = 1:size(rValues,2)      
            if(stable(iKappa,iUnique))
                colour = [0,0.7,0];
            else
                colour = [0.7,0,0];
            end
%             plot(kappaList(iKappa),rValues(iUnique,iR),...
%                 'marker',markerList{iR}, ...
%                 'markeredgecolor',colour,'markersize',10,'linewidth',2);

            plot(kappaList(iKappa),rValues(iUnique,iR),...
                'marker',markerList{iUnique}, ...
                'markeredgecolor',colour,'markersize',10,'linewidth',2);


            ylim([-1,1]);
            hold on;
            
        end
%       
        xlabel('$\kappa (\beta = 1)$', 'interpreter','latex');
        ylabel('$r_2$','interpreter','latex');


    end
end

return





