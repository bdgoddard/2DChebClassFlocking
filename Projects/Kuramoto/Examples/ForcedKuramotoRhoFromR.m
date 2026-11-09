function rho = ForcedKuramotoRhoFromR(rc,rs,aLine,opts)
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
    q = aLine.Pts.y;
    % integration vector
    Int = aLine.Int;

    % parameters from the kernel
    a = opts.a;   % a vector of cosine fourier coefficients [a_1;a_2;...;a_n]
    kappa = opts.kappa; % strength of the interaction
    beta = opts.beta;   % the diffusion term has \beta^{-1} in front of it
    
    % parameters for the potential
    h = opts.h; % strength of the potential, h \cos\theta

    % parameters for the forcing
    F = opts.F; % strength of the forcing, F \theta

    % number of Fourier modes
    n = length(a);

    % precompute cos(k*q) and sin(k*q); this isn't really necessary, but we
    % use them a lot
    C = zeros(length(q),n);
    S = zeros(length(q),n);
    for k = 1:n
        C(:,k) = cos(k*q);
        S(:,k) = sin(k*q);
    end

    % construct the exponent as the linear combination of 
    % a_k ( r_c,k \cos(k q) + r_s,k \sin(k q) )
    exponent = 0;
    for j = 1:n
        exponent = exponent - a(j) * ( rc(j)*C(:,j) + rs(j)*S(:,j) ) ;
    end
    
    % multiply by the strength
    exponent = kappa * exponent;

    % add in the potential
    exponent = exponent + h * C(:,1);

    % multiply the whole thing by beta
    exponent = beta * exponent;

    % compute the integral (which depends on q)
    integral = zeros(size(q));
    nPts = length(q);
    for iPt = 1:nPts

        % form the exponent term from the (shifted) interaction
        integralExponent = 0;
        theta = q(iPt);
        for j = 1:n
            integralExponent = integralExponent - a(j) * ( rc(j)*cos(j*(q+theta)) + rs(j)*sin(j*(q+theta)) ) ;
        end
        integralExponent = kappa * integralExponent;

        % add in the (shifted) potential
        integralExponent = integralExponent + h * cos(q+theta);

        % add in forcing, both the term from the integrand and the one
        % from the exponent to enforce periodicity
        integralExponent = integralExponent + F * q;

        % apply -\beta factor
        integralExponent = - beta * integralExponent;

        integral(iPt) = Int * exp(integralExponent);
    end

    % compute the r's by doing the weighted integration and
    % normalizing
    rho = exp(exponent) .* integral;
    Z = Int * rho;
    rho = rho / Z;

end
