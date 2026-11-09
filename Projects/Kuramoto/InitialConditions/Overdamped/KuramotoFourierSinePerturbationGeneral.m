function rho0 = KuramotoFourierSinePerturbationGeneral(q,opts)
    
    rho0 = KuramotoRhoFourier(q,opts);
    
    if(isfield(opts,'b'))
        b = opts.b;
    else
        b = 1;
    end
    
    nMax = length(b);
    
    if(isfield(opts,'shift'))
        shift = opts.shift;
    else
        shift = 0;
    end
    
    perturb = zeros(size(rho0));
    for j = 1:nMax
        perturb = perturb + b(j)*sin(j*q + shift);
    end    
    
    epsilon = opts.epsilon;
    rho0 = rho0.*(1 + epsilon*perturb);
    
end