function data = KuramotoOverdamped(opts,optsPlot)

% Ensure we have dirData as a variable
AddPaths;
% Set up save location
saveDirBase = ['Kuramoto' filesep 'Overdamped'];

if(nargin<1)
    opts = struct;
end
if(nargin<2)
    optsPlot = struct;
    optsPlot.doPlots = true;
end
    
%--------------------------------------------------------------------------
% Set some sensible defaults
%--------------------------------------------------------------------------

default = KuramotoDefaults();

%--------------------------------------------------------------------------
% Combine input options with defaults
%--------------------------------------------------------------------------

props = {'geom','geomPlot','doLog','beta','kappa','relTol','absTol', ...
         'optsIC','optsKernel','optsVext','tMax','nSaves','toEq','eqTol'};
     
nProps = length(props);
for iProp = 1:nProps
    prop = props{iProp};
    if(~isfield(opts,prop))
        opts.(prop) = default.(prop);
    end
end

% set phase space flag
opts.phasespace = false;

%--------------------------------------------------------------------------
% Set up domain and get relevant options
%--------------------------------------------------------------------------

geom = opts.geom;
geomPlot = opts.geomPlot;

aLine = FourierLine(geom);
aLine.ComputeAll(geomPlot);

L = aLine.L;

tMax = opts.tMax;
nSaves = opts.nSaves;
plotTimes = (0:tMax/nSaves:tMax)';

beta = opts.beta;
kappa = opts.kappa;
optsKernel = opts.optsKernel;
optsKernel.L = L;
optsVext = opts.optsVext;
optsVext.L = L;

relTol = opts.relTol;
absTol = opts.absTol;

toEq = opts.toEq;
if(~toEq)
    opts.eqTol = NaN;
end    
eqTol = opts.eqTol;

doLog = opts.doLog;

optsIC = opts.optsIC;
optsIC.L = L;

% parameters we care about
optsPhys = v2struct(beta,kappa,tMax,nSaves,optsIC,optsKernel,optsVext);
optsNum  = v2struct(geom,doLog,relTol,absTol,toEq,eqTol);
optsPDE = v2struct(optsNum,optsPhys);

% Modify save location to make it quicker to look up computations
saveDir = [saveDirBase filesep 'kappa_' num2str(kappa)];

%----------------------------------------------------------------------
% Initialisation
%----------------------------------------------------------------------

q = aLine.Pts.y;
Dq = aLine.Diff.Dy;
DDq = aLine.Diff.DDy;
Int = aLine.Int;

convOpts.geom = opts.geom;
convOpts.optsKernel = optsKernel; % note this is updated with L

% convOpts.optsKernel.a_n

convOpts.convDir = [saveDirBase filesep 'Convolutions' filesep convOpts.optsKernel.kernelFn];
Conv = ComputeConvolutionWithParameters(aLine,convOpts);
aLine.Conv = Conv;

% % check of an exact formula with a_n = -1 (-cos(q))
% alpha = 2*pi/L;
% convTest = Conv*(sin(alpha*q).^3) + 3/4*pi/alpha*sin(alpha*q);
% optsPlot.plain = true;
% optsPlot.linewidth = 1.5;
% optsPlot.linecolor = 'r';
% aLine.plot(Conv*(sin(alpha*q).^3),optsPlot)
% hold on
% optsPlot.linecolor = 'b';
% optsPlot.linestyle = '--';
% aLine.plot(-3/4*pi/alpha*sin(alpha*q),optsPlot)
% max(abs(convTest))
% return

if(doLog)
    rhsFn = @dlogRho_dt;
else
    rhsFn = @drho_dt;
end

mM    = ones(size(q));

rho_ic = KuramotoIC(q,optsIC,Int);

if(doLog)
    soln_ic = log(rho_ic);
else
    soln_ic = rho_ic;
end

VextFn = str2func(optsVext.VextFn);
V = VextFn(q,optsVext);

% figure
% optsPlot.plain = true;
% aLine.plot(Dq*rho_ic,optsPlot);
% hold on
% optsPlot.linestyle = '--';
% aLine.plot(2*pi/optsIC.L*optsIC.epsilon*cos(2*pi/optsIC.L*q+optsIC.shift),optsPlot);
% 
% Int*rho_ic
% pause


%---------------------------------------------------------------------%
% Solve the PDE                                                       %
%---------------------------------------------------------------------%

optsOther.SkipLoadQuestion = true;
recompute = false;

[data,~,paramsDS] = DataStorage(saveDir,...
                            @ComputeDynamics,optsPDE,optsOther,[],[],recompute);

if(optsPlot.doPlots)
    rho_t = data.rho_t;
    
    rhoMax = 1.1*max(rho_t(:));

    optsPlot.plain = true;
    optsPlot.linecolor = 'r';
    optsPlot.linewidth = 1.5;

    optsPicard.V = V;
    optsPicard.Conv = Conv;
    optsPicard.kappa = kappa;
    optsPicard.beta = beta;
    
    rhoInf = rhoKuramotoPicard(rho_t(:,end),aLine,optsPicard);
    optsPlotInf = optsPlot;
    optsPlotInf.linecolor = 'k';
    optsPlotInf.linestyle = '--';

    figure
    for iT = 1:length(plotTimes)
        rho_i = rho_t(:,iT);
        hold off
        aLine.plot(rho_i,optsPlot);
        hold on
        aLine.plot(rhoInf,optsPlotInf);
        ylim([0,rhoMax]);
        title(plotTimes(iT));
        pause(0.1)
    end
    
    m = data.m;
    disp(['Maximum relative mass error: ' num2str(max(abs(m-m(1))/m(1))) ]);
    rhoErr = Int*(abs(rhoInf - rho_t(:,end)));
    disp(['Maximum L1 error against Picard: ' num2str(rhoErr) ]);
end
    
%----------------------------------------------------------------------
% RHS of PDE
%----------------------------------------------------------------------

function drhodt = drho_dt(t,rho) 

    drhodt = beta^-1 * DDq*rho + Dq*( rho.*( Dq*(kappa*Conv*rho + V) ) );
    
end

function dlogRhodt = dlogRho_dt(t,logRho) 

    % rho = exp( h(q,t) )
    % drho/dt = dh/dt * rho
    % drho/dx = dh/dx * rho
    % d^2rho/dx^2 = d^2h/dx^2 * rho + (dh/dx)^2*rho
        
    rho = exp(logRho);
    potentialTerm = kappa*Conv*rho + V;
    
    dlogRhodt = beta^-1 * ( DDq*logRho + (Dq*logRho).^2 ) ...
                 +  ( (Dq*logRho).*(Dq*potentialTerm) ...
                             + DDq*potentialTerm );
end

%----------------------------------------------------------------------
% Wrapper for dynamics
%----------------------------------------------------------------------

function data = ComputeDynamics(params,misc)                        
        
    if(toEq)
        optsODE = odeset('RelTol',relTol,'AbsTol',absTol,'Mass',diag(mM),'Events',@eventfun);
        %[~,logRho_t,tEq] = ode15s(rhsFn,plotTimes,logRho_ic,optsODE);
        [~,soln_t,tEq] = ode15s(rhsFn,plotTimes,soln_ic,optsODE);
        if(isempty(tEq))
            tEq = NaN;
        end
    else
        optsODE = odeset('RelTol',relTol,'AbsTol',absTol,'Mass',diag(mM));
        [~,soln_t] = ode15s(rhsFn,plotTimes,soln_ic,optsODE);
        tEq = NaN;
    end
    
    soln_t = soln_t';
    
    if(doLog)
        rho_t = exp(soln_t);
    else
        rho_t = soln_t;
    end
    
    if(toEq)
        nBeforeEq = size(rho_t,2);
        nSaves = length(plotTimes);
        if(nBeforeEq < nSaves)
            rho_t_Full = zeros(size(rho_t,1),length(plotTimes));
            rho_t_Full(:,1:nBeforeEq) = rho_t;
            for iT = nBeforeEq+1:nSaves
                rho_t_Full(:,iT) = rho_t(:,nBeforeEq);
            end
            rho_t = rho_t_Full;
        end            
    else
        nSaves = length(plotTimes);
    end

    F = zeros(nSaves,1);
    optsFE.beta = opts.beta;
    optsFE.kappa = opts.kappa;
    optsFE.optsVext = opts.optsVext;
    
    m = zeros(nSaves,1);
    
    for iSave = 1:nSaves
        F(iSave) = KuramotoFreeEnergy(rho_t(:,iSave),aLine,optsFE);
        m(iSave) = Int*rho_t(:,iSave);
    end
    
    data = v2struct(opts,aLine,rho_t,plotTimes,F,m,tEq);
    
end

%----------------------------------------------------------------------
% Equilibrium event function
%----------------------------------------------------------------------

    function [err,isterm,dir] = eventfun(t,rho)
        drhodt = rhsFn(t,rho);
        %t
        err = max(abs(drhodt)) - eqTol;
        isterm = 1;
        dir = 0;
    end

end