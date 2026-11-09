function rho = KuramotoIC(q,opts,Int)
    
    ICfn = str2func(opts.rhoICFn);
    
    rho = ICfn(q,opts);
    
    if(nargin>2)
        Z = Int*rho;
        rho = rho/Z;
    end
    
end
