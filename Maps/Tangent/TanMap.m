function [z,dz,dx,ddx,dddx,ddddx,ddz] = TanMap(x,alpha1,alpha2)
%[z,dz,dx,ddx,dddx] = TanMap(x,alpha1,alpha2)
% z     = alpha2 + tan( (x -s0)*lambda )/alpha1;
% [-1,1] -> [-1,1]
% if alpha1>>1, points cluster around alpha2
% inverse of AtanMap
% A. Bayliss and E. Turkel. Mappings and accuracy for Cheby-
% shev pseudospectral approximations. Journal of Computational
% Physics, 101:349{359, 1992.

    kappa  = atan(alpha1.*(1+alpha2))./atan(alpha1.*(1-alpha2));
    s0     = (kappa - 1)./(kappa + 1);
    lambda = atan(alpha1.*(1-alpha2))./(1-s0);

    z  = alpha2 + tan( (x - s0).*lambda )./alpha1;
    X  = (x - s0).*lambda;
    dz = lambda./alpha1 .* sec(X).^2;

    ddz = 2 * lambda .* dz .* tan( (x - s0).*lambda );
    
    %Inverse
    
    %x     = s0 + atan(alpha1*(z-alpha2))/lambda;
    k      = (alpha1.^2.*(z-alpha2).^2 + 1);
    dx     = alpha1./lambda ./ k;

    ddx    = - 2*alpha1./lambda .* alpha1.^2.*(z-alpha2) ./ k.^2;
    
    dddx   = 2*alpha1.^3./lambda .* (3*alpha1.^2.*(z-alpha2).^2 - 1) ./ k.^3;
    
    ddddx  = - 24*alpha1.^5/lambda .* (z-alpha2).*(alpha1.^2.*(z-alpha2).^2 - 1) ./ k.^4;

end