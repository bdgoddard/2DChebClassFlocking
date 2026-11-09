function rho = KuramotoRhoFourier(q,opts)

    a = opts.a;
    r = opts.r;
    beta = opts.beta;
    kappa = opts.kappa;
    
    n = length(a);

    exponent = 0;

    for j = 1:n
        exponent = exponent + a(j)*r(j)*cos(j*q);
    end
    exponent = -beta*kappa*exponent;

    rho = exp(exponent);

end
