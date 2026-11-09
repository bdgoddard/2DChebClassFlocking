function [aLine,Rho_t] = ForcedKuramotoPDETest(opts)

if(nargin == 0)
    opts = struct;
end

% choose the number of points
if(isfield(opts,'N'))
    geom.N = opts.N;
else
    geom.N = 30;
end

% domain is [0,2pi]
geom.yMin = 0;
geom.yMax = 2*pi;

% set parameters
if(isfield(opts,'D0')) % diffusion parameter
    D0 = opts.D0;
else
    D0 = 1;
end

if(isfield(opts,'F')) % forcing parameter
    F = opts.F;
else
    F = 1;
end

if(isfield(opts,'kappa')) % interaction parameter
    kappa = opts.kappa;
else
    kappa = 2.5;
end

% maximum time
% !! change this depending on how close you are to the critical kappa
if(isfield(opts,'tMax'))
    tMax = opts.tMax;
else
    tMax = 10;
end

% set up the domain
aLine = FourierLine(geom);

% compute differentiation and integration
[Pts,Diff,Int,~] = aLine.ComputeAll();
% interpolation for plotting
Interp             = aLine.ComputeInterpolationMatrix((0:0.01:1)',true);

% compute the convolution term
Conv = aLine.ComputeConvolutionMatrix(@KernelWrapper,true);

% convOpts.optsKernel.kappa = kappa;
%
% Conv = ComputeConvolutionWithParameters(aLine,convOpts);

% shorthand the differentiation matrix
Dy = Diff.Dy;

% set up the initial condition
rho_ic = InitialCondition(Pts.y,0);

% times at which to output the results
outTimes = [0:tMax/100:tMax];

% solve the PDE (periodic boundary conditions are done automatically)
opts = odeset('RelTol',10^-9,'AbsTol',10^-9);
[~,Rho_t] = ode15s(@rhs,outTimes,rho_ic,opts);

% plot the results
maxRho = max(max(Rho_t));
figure
for iTimes = 1:length(outTimes)
    rhot = Rho_t(iTimes,:)';

    hold off
    aLine.plot(rhot);

    ylim([0,maxRho]);

    pause(0.05)

end

%----------------------------------------------------------------------

    function drhodt = rhs(t,rho)

        % get the flux
        flux = getFlux(rho);

        % set up the PDE
        drhodt = -Dy * flux;

        drhodt = drhodt(:);

    end

    function f = getFlux(rho)
        % get the flux from teh three terms
        j = getAdvection(rho);
        f2 = get2Body(rho);

        % note that the flux has a minus sign by convention
        f = - D0*Dy*rho - j - f2;
    end

    function j = getAdvection(rho)
        % compute the advection term
        j = - F * rho;
    end

    function f = get2Body(rho)
        % compute the convolution term
        f = rho.*(Dy*(Conv * rho));
    end


    function y = InitialCondition(y,t)
        % essentially arbitrarily-chosen periodic intitial condition
        y = (sin(y) + 2);
        y = y/(Int*y);  % ensure it's normalised

    end


% set up the kernel in a slightly odd way so it works with the
% convolution code
    function g = KernelWrapper(y)
        g = Kernel(y,kappa);
    end

    function g = Kernel(y,kappa)
        g =  - kappa * cos(y);
    end

end