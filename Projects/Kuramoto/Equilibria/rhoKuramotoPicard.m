function rho = rhoKuramotoPicard(rhoCurrent,aLine,opts)

    % these parameters are largely made up
    if(isfield(opts,'alpha'))
        alpha = opts.alpha;
    else
        alpha = 0.1;
    end
    if(isfield(opts,'tol'))
        tol = opts.Tol;
    else
        tol = 10^(-7);
    end

    if(~isempty(aLine))
        Int = aLine.Int;
    else
        Int = opts.Int;
    end
       
    V = opts.V;
    Conv = opts.Conv;
    
    beta = opts.beta;
    kappa = opts.kappa;
    
    % intitial error is always higher than the tolerance
    err = Inf;

    Z = Int*rhoCurrent;
    rhoCurrent = rhoCurrent/Z;

    count = 0;

    while (err > tol)

        count = count + 1;

        rhoOld = rhoCurrent;

        rhoNew = exp(- beta*( V + kappa*Conv*rhoCurrent ));
        Z = Int*rhoNew;
        rhoNew = rhoNew/Z;

        rhoCurrent = (1-alpha)*rhoCurrent + alpha*rhoNew;

        err = Int*(abs(rhoCurrent - rhoOld));

    end

    rho = rhoCurrent;

end
