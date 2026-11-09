function KuramotoCriticalTemperature(optsIn,optsPert)

%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------

AddPaths;
% Set up save location
saveDirBase = ['Kuramoto' filesep 'Overdamped' filesep 'CriticalTemperature'];

%--------------------------------------------------------------------------
% Load some sensible defaults
%--------------------------------------------------------------------------

opts = KuramotoDefaults();

if(isfield(optsIn,'geom'))
    opts.geom = optsIn.geom;
end

%--------------------------------------------------------------------------
% Kernel options
%--------------------------------------------------------------------------
%
% Note that this should work for a general kernel, of arbitrary length.  We
% are setting the fourier modes of \sum a_n cos(nq) in a vector
% [a_1; a_2; ...; a_n]

opts.kappa = 1; % we're going to vary beta
opts.optsKernel.a_n = optsIn.optsKernel.a_n;
opts.beta = optsIn.beta;

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
% Compute stability
%--------------------------------------------------------------------------

fprintf('Finding solutions to r equations ...\n');

% use DataStorage for rUnique
recompute = false;
optsOther.SkipLoadQuestion = true;
optsOther.aLine = aLine;
rUniqueDir = [saveDirBase filesep 'r'];

optsRUnique = opts;
optsRUnique.dr = 0.1; % so we can control this

dataRUnique = DataStorage(rUniqueDir,@getRUnique,optsRUnique,optsOther,[],[],recompute);

rUnique = dataRUnique.rUnique;

fprintf('Testing stability ...\n');
testStability(optsIn,optsPert,rUnique)


%-------------------------------------------------------------------------- 
    
    function data = getRUnique(opts,optsOther)
        %--------------------------------------------------------------------------
        % Set up a grid of initial guesses for r, which we hope will lead to all
        % possible solutions of the self-consistency equations
        %--------------------------------------------------------------------------
        
        aLine = optsOther.aLine;
        
        % the 1D grid for each r (we could generalize this to different values, but
        % this seems good enough)        
        rGrid = -1:opts.dr:1;

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
            disp(iTest/nTests);
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
        [rUnique,rIndex] = uniquetol(rTests,uTol,'ByRows',true);

        % display the results
        rUnique = sortrows(rUnique,'descend');
        
        data.rUnique = rUnique;
    end
    
    function testStability(optsIn,optsPert,rUnique)

        %--------------------------------------------------------------------------
        % Perturb the various fixed point rho's and run the PDE to check for
        % stability
        %--------------------------------------------------------------------------

        opts = optsIn;
        
        % overwrite initial condition to the perturbed one
        opts.optsIC = optsPert.optsIC;
       
        % set the final time - this is highly dependent on the choice of parameters        
        opts.tMax = optsPert.tMax;

        % set the number of save points (to make a nice movie and plots of the free energy)
        opts.nSaves = 50;

        % turn off the plotting
        optsPlot.doPlots = false;

        if(optsPlot.doPlots)
            figure
        end

        % run the PDE dynamics using each value of r
        nUnique = size(rUnique,1);
        for iUnique = 1:nUnique
            opts.optsIC.r = rUnique(iUnique,:)';
            data(iUnique) = KuramotoOverdamped(opts,optsPlot);
            pause(1) % ensure different file names
        end

        % plot the free energies in the same colours as the rho plots above
        hf = figure('Position',[100,100,1200,800]);
        maxF = -inf;
        minF = inf;

        nDim = length(opts.optsKernel.a_n);
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
    end

    function plotRhos(opts,rUnique)
    
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
%             test_r_i = rKuramotoPicard(r_i,aLine,optsPicard)';
            hold on    
        end    
        legend
    end




end