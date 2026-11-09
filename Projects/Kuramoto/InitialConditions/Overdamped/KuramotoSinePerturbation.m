function rho0 = KuramotoSinePerturbation(q,opts)

    if(~isfield(opts,'L'))
        L = 2*pi;
    else
        L = opts.L;
    end
    alpha = 2*pi/L;

    epsilon = opts.epsilon;
    shift = opts.shift;
    rho0 = ones(size(q)) + epsilon*sin(alpha*q+shift);
end