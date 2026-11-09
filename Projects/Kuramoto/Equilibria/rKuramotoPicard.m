function r = rKuramotoPicard(rCurrent,aLine,opts)
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

    % precompute cos(k*q); this isn't really necessary
    C = zeros(length(q),n);
    for k = 1:n
        C(:,k) = cos(k*q);
    end

    % preallocate the solution
    rNew = zeros(n,1);

    % perform a Picard iteration until the error is less than the set
    % tolerance
    while (err > tol)

        % keep track of the old r's to compute the error with
        rOld = rCurrent;

        % form the new r's, which are given by the right hand side of the
        % self-consistency equations
        for k = 1:n
            % construct the exponent as the linear combination of 
            % a_k r_k \cos(k q)
            exponent = 0;
            for j = 1:n
                exponent = exponent - a(j)*rCurrent(j)*C(:,j);
            end
            
            % note that it's actually only the product of \kappa and \beta
            % that matters; you can see this by multiplying the PDE by
            % \beta, for example
            exponent = beta*kappa*exponent;

            % compute the r's by doing the weighted integration and
            % normalizing
            rNew(k) = Int*(C(:,k).*exp(exponent)) ...
                      / ( Int*(exp(exponent)) );
        end

        % apply the (weighted) Picard scheme
        rCurrent = (1-alpha)*rCurrent + alpha*rNew;
        
        % compute the error to see if we've converged
        err = max(abs(rCurrent - rOld));

    end
    
    % return the computed values
    r = rCurrent;

end
