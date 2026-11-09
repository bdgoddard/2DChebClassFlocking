function f0 = KuramotoFourierSinePerturbationMB(q,p,qInt,pInt,opts)
            
    phi0 = KuramotoMB(p,pInt,opts);
    
    rho0 = KuramotoFourierSinePerturbation(q,qInt,opts);
    
    f0 = rho0.*phi0;
    
end