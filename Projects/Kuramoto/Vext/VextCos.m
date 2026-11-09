function [V,DV] = VextCos(q,opts)
    alpha = opts.alpha;
    phi = opts.phi;
    
    V = -alpha*cos(q+phi);
    
    if(nargout == 2)
        DV = alpha*sin(q+phi);
    end
end