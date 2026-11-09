function [rc,rs] = rForcedKuramotoPicard(rcCurrent,rsCurrent,aLine,opts)
% INPUTS: 
% - rCurrent is a vector of the initial guesses of the n r_j to start
% the Picard scheme.
% - aLine contains points and integration weights; you could just pass
% these separately
% - opts contains the a_j, \kappa, and \beta paramters, and optionally the
% mixing paramter alpha for the Picard scheme

    % This uses the pseudospectral machinery from 2DChebClass.  However,
    % this isn't really necessary.  You just need a grid of points and a
    % way to compute integrals reasonably accurately - finite difference is
    % probably fine.
    
    % points (this is theta for the Kuramoto model)
    q    = aLine.Pts.y;
    nPts = aLine.N;
    % integration vector
    Int  = aLine.Int;

     % parameters from the kernel
    a = opts.a;   % a vector of cosine fourier coefficients [a_1;a_2;...;a_n]
    kappa = opts.kappa; % strength of the interaction
    beta = opts.beta;   % the diffusion term has \beta^{-1} in front of it

    % parameters for the potential
    h = opts.h; % strength of the potential, h \cos\theta

    % parameters for the forcing
    F = opts.F; % strength of the forcing, F \theta

    % tolerance to solve to
    tol = opts.tol;

 % a mixing parameter - this stabilises the Picard iteration scheme
    if(isfield(opts,'alpha'))
        alpha = opts.alpha;
    else
        % this parameter is largely made up
        alpha = 0.1;
    end

    % intitial error is always higher than the tolerance
    err = Inf;

    % number of Fourier modes
    n = length(a);

    % -------- Phase-fixing options ----------
    % Default: enable phase fixing only when |h| is tiny
    if isfield(opts,'phaseFix')
        phaseFix = opts.phaseFix;
    else
        phaseFix = (abs(h) < 1e-14);
    end

    % Never phase-fix when h is not ~0 (avoid changing behaviour)
    if abs(h) >= 1e-14
        phaseFix = false;
    end

    if isfield(opts,'mu_target')
        mu_target = opts.mu_target;
    else
        mu_target = pi;   % sets where you want to fix the mean
    end

    if isfield(opts,'phaseTol')
        phaseTol = opts.phaseTol;
    else
        phaseTol = 1e-8;  % robust default; increase if mu gets noisy near transition
    end
    % --------------------------------------------------------------------

    % precompute cos(k*q) and sin(k*q); this isn't really necessary, but we
    % use them a lot
    C = zeros(length(q),n);
    S = zeros(length(q),n);
    for k = 1:n
        C(:,k) = cos(k*q);
        S(:,k) = sin(k*q);
    end

    % preallocate the solution
    rcNew = zeros(n,1);
    rsNew = zeros(n,1);

    while (err > tol)

        % Store "old" values in the SAME gauge that we will use to measure error.
        % If phase-fixing is on, we keep everything gauge-fixed at each iteration.
        rcOld = rcCurrent;
        rsOld = rsCurrent;

        % form the new r's, which are given by the right hand side of the
        % self-consistency equations
        % loop over the Fourier modes
        
        % construct the exponent as the linear combination of 
        % a_k ( r_c,k \cos(k q) + r_s,k \sin(k q) )

        exponent = 0;
        for j = 1:n
            exponent = exponent - a(j) * ( rcCurrent(j)*C(:,j) + rsCurrent(j)*S(:,j) );
        end

        
        % multiply by the strength
        exponent = kappa * exponent;

        % add in the potential
        exponent = exponent + h * C(:,1);

        % multiply the whole thing by beta
        exponent = beta * exponent;

        % Compute the integral (which depends on q)
        integral = zeros(size(q));
        for iPt = 1:nPts
            
            % form the exponent term from the (shifted) interaction
            integralExponent = 0;
            theta = q(iPt);

            for j = 1:n
                integralExponent = integralExponent ...
                    - a(j) * ( rcCurrent(j)*cos(j*(q+theta)) + rsCurrent(j)*sin(j*(q+theta)) );
            end

            integralExponent = kappa * integralExponent;
            % add in the (shifted) potential
            integralExponent = integralExponent + h * cos(q+theta);

            % add in forcing, both the term from the integrand and the one
            % from the exponent to enforce periodicity
            integralExponent = integralExponent + F * q;

            % apply -beta factor
            integralExponent = -beta * integralExponent;

            integral(iPt) = Int * exp(integralExponent);
        end

        % compute the r's by doing the weighted integration and
        % normalizing
        rho = exp(exponent) .* integral;
        Z = Int * rho;
        rho = rho / Z;

        for k = 1:n
            rcNew(k) = Int*(C(:,k) .* rho);
            rsNew(k) = Int*(S(:,k) .* rho);
        end

        
        % apply the (weighted) Picard scheme
        rcCurrent = (1-alpha)*rcCurrent + alpha*rcNew;
        rsCurrent = (1-alpha)*rsCurrent + alpha*rsNew;

        % ---------------- Phase fixing (ONLY when h ~ 0) ----------------
        % Pin mean direction mu = arg(m1) to mu_target.
        if phaseFix
            m1mag = hypot(rcCurrent(1), rsCurrent(1)); % we check whether r_c^2 + r_s^2 is far enough from 0
            if m1mag > phaseTol
                mu = atan2(rsCurrent(1), rcCurrent(1)); % atan2(y,x) - this is the formula for the mean that you can find analytically
                delta = mu_target - mu;

                % Rotate every Fourier mode by k*delta
                for k = 1:n
                    ck = cos(k*delta);
                    sk = sin(k*delta);
                    rc_k = rcCurrent(k);
                    rs_k = rsCurrent(k);
                    rcCurrent(k) = rc_k*ck - rs_k*sk;
                    rsCurrent(k) = rc_k*sk + rs_k*ck;
                end

                % Also rotate the stored "old" values by the same delta,
                % so the error is measured in a consistent gauge.
                for k = 1:n
                    ck = cos(k*delta);
                    sk = sin(k*delta);
                    rc_k = rcOld(k);
                    rs_k = rsOld(k);
                    rcOld(k) = rc_k*ck - rs_k*sk;
                    rsOld(k) = rc_k*sk + rs_k*ck;
                end
            end
        end
        % ----------------------------------------------------------------

        % ---------------- Convergence check ----------------
        errc = max(abs(rcCurrent - rcOld));
        errs = max(abs(rsCurrent - rsOld));
        err  = max(errc, errs);

    end

    rc = rcCurrent;
    rs = rsCurrent;

end