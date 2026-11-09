function [K,DK] = KuramotoKernelFourier(y,opts)
% returns the kernel and its derivative (for stochastics)

    if(~isfield(opts,'L'))
        L = 2*pi;
    else
        L = opts.L;
    end
    alpha = 2*pi/L;

    a_n = opts.a_n;
    nModes = length(a_n);
    K = zeros(size(y));
    DK = zeros(size(y));
    for n = 1:nModes
        K = K + a_n(n).*cos(alpha*n*y);
        DK = DK - alpha*n*a_n(n).*sin(alpha*n*y);
    end
end
