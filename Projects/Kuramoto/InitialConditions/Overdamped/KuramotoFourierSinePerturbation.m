function rho0 = KuramotoFourierSinePerturbation(q,opts)
    
    rho0 = KuramotoRhoFourier(q,opts);
    
    if(isfield(opts,'n'))
        n = opts.n;
    else
        n = 1;
    end
    
    epsilon = opts.epsilon;
    shift = opts.shift;
    rho0 = rho0.*(1 + epsilon*sin(n*q+shift));
    
end