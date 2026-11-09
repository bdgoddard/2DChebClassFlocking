function [K,DK] = KuramotoKernelCosSquared(y,opts)
% returns the kernel and its derivative (for stochastics)

    a_n = opts.a_n;
    
    K = a_n * cos(y).^2;
    DK = -2 * a_n * cos(y).*sin(y);

end
