function f0 = KuramotoSinePerturbationMB(q,p,Int,opts)
            
    % not actually worried about the integration here as we re-normalize
    % later anyway
    phi0 = KuramotoMB(p,Int,opts);
    rho0 = KuramotoSinePerturbation(q,Int,opts);
        
    f0 = rho0.*phi0;
    Z = Int*f0;
    f0 = f0/Z;
    
end