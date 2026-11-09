function data = KuramotoPhaseSpace(opts,optsPlot)

% need to generalize to Fourier potentials with other steady states, etc

% Ensure we have dirData as a variable
AddPaths;
global dirData
% Set up save location
saveDirBase = [dirData filesep 'Kuramoto' filesep 'PhaseSpace'];

if(nargin<1)
    opts = struct;
end
if(nargin<2)
    optsPlot = struct;
    optsPlot.doPlots = true;
end

%--------------------------------------------------------------------------
% Combine input options with defaults
%--------------------------------------------------------------------------

default = KuramotoDefaultsPhaseSpace();

%--------------------------------------------------------------------------
% Combine input options with defaults
%--------------------------------------------------------------------------

props = {'geom','geomPlot','beta','kappa','gamma','relTol','absTol', ...
         'optsIC','optsKernel','optsVext','tMax','nSaves','toEq','eqTol'};

nProps = length(props);
for iProp = 1:nProps
    prop = props{iProp};
    if(~isfield(opts,prop))
        opts.(prop) = default.(prop);
    end
end

% copy beta for MB distribution
opts.optsIC.beta = opts.beta;

% set phase space flag
opts.phasespace = true;

%--------------------------------------------------------------------------
% Set up domain and get relevant options
%--------------------------------------------------------------------------

geom = opts.geom;
geomPlot = opts.geomPlot;

aBox = PeriodicBox(opts.geom);
aBox.ComputeAll(opts.geomPlot);

N1 = aBox.N1;
N2 = aBox.N2;
onesN1 = ones(N1,1);
onesN2 = ones(N2,1);

tMax = opts.tMax;
nSaves = opts.nSaves;
plotTimes = (0:tMax/nSaves:tMax)';

gamma = opts.gamma;
beta = opts.beta;
kappa = opts.kappa;

optsKernel = opts.optsKernel;
optsVext = opts.optsVext;

relTol = opts.relTol;
absTol = opts.absTol;

toEq = opts.toEq;
if(~toEq)
    opts.eqTol = NaN;
end    
eqTol = opts.eqTol;

optsIC = opts.optsIC;

% parameters we care about
optsPhys = v2struct(beta,kappa,tMax,nSaves,optsIC,optsKernel,optsVext);
optsNum  = v2struct(geom,relTol,absTol,toEq,eqTol);
optsPDE = v2struct(optsNum,optsPhys);

% Modify save location to make it quicker to look up computations
saveDir = [saveDirBase filesep 'kappa_' num2str(kappa) '_gamma_' num2str(gamma)];

%--------------------------------------------------------------------------
% Initialisation
%--------------------------------------------------------------------------


q = aBox.Pts.y2_kv;
p = aBox.Pts.y1_kv;
Dq = aBox.Diff.Dy2;
Dp = aBox.Diff.Dy1;
Int = aBox.Int;

% 1D integration vectors for computing traces
[pInt,qInt] = aBox.Compute1DIntegrationMatrices;

% 1D convolution to apply only to rho
convOpts.geom = opts.geom;
convOpts.optsKernel = opts.optsKernel;
convOpts.convDir = [saveDirBase filesep 'Convolutions' filesep convOpts.optsKernel.kernelFn];
convOpts.convFn = 'ComputeConvolutionMatrix_1D_2';
Conv_1D = ComputeConvolutionWithParameters(aBox,convOpts);

% Set up ODE function and mass matrix
rhsFn = @df_dt;
mM         = ones(size(q));
% boundaries in p (periodic in q)
pMinBound = aBox.Ind.left;
pMaxBound = aBox.Ind.right;
mM(pMinBound) = 0;
mM(pMaxBound) = 0;

% Initial condition
f_ic = KuramotoPhaseSpaceIC(q,p,optsIC,Int);

% external potential
VextFn = str2func(optsVext.VextFn);
V = VextFn(q,optsVext);

%--------------------------------------------------------------------------
% Solve the PDE                                                       
%--------------------------------------------------------------------------

optsOther.SkipLoadQuestion = true;
recompute = false;

[data,~,paramsDS] = DataStorage(saveDir,...
                         @ComputeDynamics,optsPDE,optsOther,[],[],recompute);
                     
%--------------------------------------------------------------------------
% Plotting                                                       
%--------------------------------------------------------------------------


if(optsPlot.doPlots)

    f_t = data.f_t;
    rho_t = data.rho_t;
    phi_t = data.phi_t;
    
    qLine = data.qLine;
    pLine = data.pLine;

    [rhoInf,phiInf] = getEquilibria();

    optsPlotInf.plain = true;
    optsPlotInf.linecolor = 'k';
    optsPlotInf.linestyle = '--';
    optsPlotInf.linewidth = 1.5;

    optsPlot.plain = true;
    optsPlot.linecolor = 'r';
    optsPlot.linewidth = 1.5;

    figure('Position',[100,100,1200,800]);

    rhoMax = 1.1*max(rho_t(:));
    phiMax = 1.1*max(phi_t(:));

    for iT = 1:length(plotTimes)
        f_i = f_t(:,iT);
        rho_i = rho_t(:,iT);
        phi_i = phi_t(:,iT);

        subplot(1,3,1)
        aBox.plot(f_i);
        xlabel('p');
        ylabel('q');

        subplot(1,3,2)
        hold off
        qLine.plot(rho_i,optsPlot);
        hold on
        qLine.plot(rhoInf,optsPlotInf);
        xlabel('q');
        ylim([0,rhoMax]);

        subplot(1,3,3)
        hold off
        pLine.plot(phi_i,optsPlot);
        hold on
        pLine.plot(phiInf,optsPlotInf);
        xlabel('p');
        ylim([0,phiMax]);

        pause(0.1)
%         pause

    end

    m = data.m;
    disp(['Maximum relative mass error: ' num2str(max(abs(m-m(1))/m(1))) ]);
    rhoErr = qLine.Int*(abs(rhoInf - rho_t(:,end)));
    disp(['Maximum L1 error against Picard: ' num2str(rhoErr) ]);
end

%----------------------------------------------------------------------
% RHS of PDE
%----------------------------------------------------------------------

    function dfdt = df_dt(t,f) 

        [rho,~] = getRhoPhi(f);
        Crho = Conv_1D*rho; % one point for each p
        Crho = replicateOverQ(Crho); % make a copy for each q

        dfdt = -Dq*(p.*f) + (Dq*(kappa*Crho + V)).*(Dp*f) ...
                + gamma*Dp*(p.*f + beta^(-1)*Dp*f);

        % Dirichlet zero in p, which is probably good enough
        dfdt(pMinBound) = f(pMinBound);
        dfdt(pMaxBound) = f(pMaxBound);

    end

    function C = replicateOverQ(Cp)
        C = sparse(kron(onesN1,Cp));
        C = C(:);
    end

%----------------------------------------------------------------------
% Wrapper for dynamics
%----------------------------------------------------------------------


    function data = ComputeDynamics(params,misc)                        

        if(toEq)
            optsODE = odeset('RelTol',relTol,'AbsTol',absTol,'Mass',diag(mM),'Events',@eventfun);
            [~,f_t,tEq] = ode15s(rhsFn,plotTimes,f_ic,optsODE);
            if(isempty(tEq))
                tEq = NaN;
            end
        else
            optsODE = odeset('RelTol',relTol,'AbsTol',absTol,'Mass',diag(mM));
            [~,f_t] = ode15s(rhsFn,plotTimes,f_ic,optsODE);
            tEq = NaN;
        end
        
        f_t = f_t.';
        
        if(toEq)
            nBeforeEq = size(f_t,2);
            nSaves = length(plotTimes);
            if(nBeforeEq < nSaves)
                f_t_Full = zeros(size(f_t,1),length(plotTimes));
                f_t_Full(:,1:nBeforeEq) = f_t;
                for iT = nBeforeEq+1:nSaves
                    f_t_Full(:,iT) = f_t(:,nBeforeEq);
                end
                f_t = f_t_Full;
            end            
        end

        m = ( Int*f_t )';
        [rho_t,phi_t] = getRhoPhi(f_t);
        
        [qLine,pLine] = getLines();
        
        data = v2struct(opts,aBox,f_t,rho_t,phi_t,m,plotTimes,qLine,pLine,tEq);

    end

%----------------------------------------------------------------------
% Equilibrium event function
%----------------------------------------------------------------------

    function [err,isterm,dir] = eventfun(t,f)
        dfdt = rhsFn(t,f);
        err = max(abs(dfdt)) - eqTol;
        isterm = 1;
        dir = 0;
    end

%----------------------------------------------------------------------
% 1D variables
%----------------------------------------------------------------------

    function [rho,phi] = getRhoPhi(f)
        rho = pInt*f;
        phi = qInt*f;
    end

%----------------------------------------------------------------------
% 1D lines
%----------------------------------------------------------------------

    function [rhoLine,phiLine] = getLines()
        
        rhoLineGeom.N = geom.N(2);
        rhoLineGeom.yMin = geom.y2Min;
        rhoLineGeom.yMax = geom.y2Max;
        rhoLine = FourierLine(rhoLineGeom);
        rhoLineGeomPlot = rhoLineGeom;
        rhoLineGeomPlot.N = 200;
        rhoLine.ComputeAll(rhoLineGeomPlot);

        phiLineGeom.N = geom.N(1);
        phiLineGeom.yMin = geom.y1Min;
        phiLineGeom.yMax = geom.y1Max;
        phiLine = SpectralLine(phiLineGeom);
        phiLineGeomPlot = phiLineGeom;
        phiLineGeomPlot.N = 200;
        phiLine.ComputeAll(phiLineGeomPlot);

    end

%----------------------------------------------------------------------
% Equilibrium traces and f
%----------------------------------------------------------------------

    function [rhoInf,phiInf] = getEquilibria()

        q = qLine.Pts.y;
        optsPicard.Conv = Conv_1D;
        optsPicard.V = VextFn(q,optsVext);
        optsPicard.kappa = kappa;
        optsPicard.beta = beta;
        rhoInf = rhoKuramotoPicard(rho_t(:,end),qLine,optsPicard);
        
        p = pLine.Pts.y;
        pInt = pLine.Int;
        phiInf = KuramotoMB(p,optsPicard);
        Z = pInt*phiInf;
        phiInf = phiInf/Z;
                        
    end

end
