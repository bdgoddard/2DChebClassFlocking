function phi0 = KuramotoMB(p,opts)
    beta = opts.beta;
    phi0 = exp(-beta*p.^2/2); 
end
