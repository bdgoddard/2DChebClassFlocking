%% sweepCriticalGammaF.m
% Parameter sweep to test whether critical gamma is independent of F
% (gamma here corresponds to 'kappa' in ForcedKuramotoPDETest)

clear; clc;


% --- Sweep settings ---
N     = 60;     % increase resolution vs 30 for clearer peaks
D0    = 1;
tMax  = 300;     % longer helps approach steady state near criticality
nOut  = 200;    % output times for solver

Fs     = [0, 0.5, 1.0];
gammas = [1.5, 2.0, 2.5, 3.5];  % example values around gamma_c ~ 2

% storage
results = struct();
results.Fs = Fs;
results.gammas = gammas;

% Solve all cases
fprintf("Running sweep...\n");
for iF = 1:numel(Fs)
    F = Fs(iF);

    for ig = 1:numel(gammas)
        gamma = gammas(ig);

        opts = struct();
        opts.N = N;
        opts.D0 = D0;
        opts.F = F;
        opts.kappa = gamma;   % <-- map gamma -> kappa used by the code
        opts.tMax = tMax;
        opts.nOut = nOut;

        [y, rhoFinal, metrics] = solveForcedKuramotoFinal(opts);

        results.rhoFinal{iF,ig} = rhoFinal;
        results.metrics(iF,ig)  = metrics;
        results.y = y;  % same grid each time (for same N)
        
        fprintf("  F = %.2f, gamma = %.2f: peak=%.4f, L2dev=%.4e\n", ...
            F, gamma, metrics.peak, metrics.L2dev);
    end
end

% --- Plot 1: final profiles for each F (each line = a gamma) ---
for iF = 1:numel(Fs)
    figure('Name', sprintf('Final profiles (F = %.2f)', Fs(iF)));
    hold on; grid on;
    for ig = 1:numel(gammas)
        plot(results.y, results.rhoFinal{iF,ig}, 'LineWidth', 1.5);
    end
    yline(1/(2*pi), '--', 'Uniform 1/(2\pi)');
    xlabel('u');
    xlim([0,6.17])
    ylabel('\rho(u,t_{max})');
    title(sprintf('Final density profiles at t_{max} (F = %.2f)', Fs(iF)));
    legend(arrayfun(@(g) sprintf('\\theta = %.2f', g), gammas, 'UniformOutput', false), ...
        'Location', 'best');
end

% --- Plot 2: deviation from uniform vs gamma for each F ---
figure('Name','Deviation from uniform vs gamma');
hold on; grid on;

for iF = 1:numel(Fs)
    L2dev = arrayfun(@(ig) results.metrics(iF,ig).L2dev, 1:numel(gammas));
    plot(gammas, L2dev, '-o', 'LineWidth', 1.5);
end
xlabel('\gamma');
ylabel('||\rho - 1/(2\pi)||_2 (discrete)');
title('Peakedness metric vs \gamma for different F');
legend(arrayfun(@(F) sprintf('F = %.2f', F), Fs, 'UniformOutput', false), 'Location', 'best');

% Optional: also show peak height above uniform
figure('Name','Peak height above uniform vs gamma');
hold on; grid on;
for iF = 1:numel(Fs)
    peakAbove = arrayfun(@(ig) results.metrics(iF,ig).peak - 1/(2*pi), 1:numel(gammas));
    plot(gammas, peakAbove, '-o', 'LineWidth', 1.5);
end
xlabel('\gamma');
ylabel('max(\rho) - 1/(2\pi)');
title('Peak height above uniform vs \gamma for different F');
legend(arrayfun(@(F) sprintf('F = %.2f', F), Fs, 'UniformOutput', false), 'Location', 'best');


%% ===== Local function: solve PDE to tMax and return final profile =====
function [y, rhoFinal, metrics] = solveForcedKuramotoFinal(opts)

    % defaults (matches your original code style)
    if ~isfield(opts,'N'),    opts.N = 30; end
    if ~isfield(opts,'D0'),   opts.D0 = 1; end
    if ~isfield(opts,'F'),    opts.F = 1; end
    if ~isfield(opts,'kappa'),opts.kappa = 1.9; end
    if ~isfield(opts,'tMax'), opts.tMax = 10; end
    if ~isfield(opts,'nOut'), opts.nOut = 100; end

    N     = opts.N;
    D0    = opts.D0;
    F     = opts.F;
    kappa = opts.kappa;
    tMax  = opts.tMax;
    nOut  = opts.nOut;

    % domain [0, 2pi]
    geom = struct();
    geom.N    = N;
    geom.yMin = 0;
    geom.yMax = 2*pi;

    % set up the domain + operators (same as ForcedKuramotoPDETest)
    aLine = FourierLine(geom);
    [Pts,Diff,Int,~] = aLine.ComputeAll();
    Dy = Diff.Dy;

    % convolution with kernel -kappa*cos(y)
    Conv = aLine.ComputeConvolutionMatrix(@KernelWrapper,true);

    % initial condition (same as your code)
    rho_ic = (sin(Pts.y) + 2);
    rho_ic = rho_ic/(Int*rho_ic); % normalize

    % output times
    outTimes = linspace(0, tMax, nOut+1);

    % solve
    odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [~,Rho_t] = ode15s(@rhs, outTimes, rho_ic, odeOpts);

    % final state
    rhoFinal = Rho_t(end,:).';
    y = Pts.y;

    % metrics for “peakedness”
    rhoUniform = 1/(2*pi);
    metrics.peak = max(rhoFinal);
    metrics.L2dev = norm(rhoFinal - rhoUniform) / sqrt(numel(rhoFinal)); % discrete L2-ish

    % ---- nested functions ----
    function drhodt = rhs(~, rho)
        flux = getFlux(rho);
        drhodt = -Dy * flux;
        drhodt = drhodt(:);
    end

    function f = getFlux(rho)
        j  = -F * rho;                 % advection part (same as original)
        f2 = rho .* (Dy*(Conv*rho));   % interaction part
        f  = -D0*Dy*rho - j - f2;      % flux (with original sign convention)
    end

    function g = KernelWrapper(yin)
        g = -kappa*cos(yin);
    end
end