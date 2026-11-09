function [V,DV] = VextCosSeries(q,opts)

    if(~isfield(opts,'L'))
        L = 2*pi;
    else
        L = opts.L;
    end
    alpha = 2*pi/L;

    % Fourier coefficients
    a_n = opts.a_n;
    
    % phase shift
    if(isfield(opts,'phi'))
        phi = opts.phi;
    else
        phi = 0;
    end
    
    % number of modes
    nModes = length(a_n);
    
    % create potential and derivative
    V = zeros(size(q));
    DV = zeros(size(q));
    
    for n = 1:nModes
        V = V + a_n(n).*cos(alpha*n*q + phi);
        DV = DV - alpha*n*a_n(n).*sin(alpha*n*q + phi);
    end

end