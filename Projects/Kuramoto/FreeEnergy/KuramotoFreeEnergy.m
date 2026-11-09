function F = KuramotoFreeEnergy(rho,aLine,opts)

    Int = aLine.Int;
    Conv = aLine.Conv;
    
    beta = opts.beta;
    kappa = opts.kappa;
    
    VextFn = str2func(opts.optsVext.VextFn);
    q = aLine.Pts.y;
    V = VextFn(q,opts.optsVext);

    F = beta^(-1) * Int*(rho .* (log(rho) - 1)) ...
        + 0.5 * kappa * Int*( rho.* (Conv*rho) ) ...
        + Int*(V.*rho);

end