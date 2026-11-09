function f = KuramotoPhaseSpaceIC(q,p,opts,Int)

    rhoICFn = str2func(opts.rhoICFn);
    phiICFn = str2func(opts.phiICFn);

    rho = rhoICFn(q,opts);
    phi = phiICFn(p,opts);
    
    f = rho.*phi;
    
    if(nargin>2)
        Z = Int*f;
        f = f/Z;
    end
    
end