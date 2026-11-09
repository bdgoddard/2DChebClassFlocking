function [z,dz,dx,ddx,dddx,ddddx] = AtanMap(x,alpha1,alpha2)
%[z,dz,dx,ddx,dddx] = AtanMap(x,alpha1,alpha2)
% z = z0 + atan(alpha1*(x-alpha2))/lambda;
% [-1,1] -> [-1,1]
% inverse of TanMap
% A. Bayliss and E. Turkel. Mappings and accuracy for Cheby-
% shev pseudospectral approximations. Journal of Computational
% Physics, 101:349{359, 1992.


    kappa  = atan(alpha1*(1+alpha2))/atan(alpha1*(1-alpha2));
    s0     = (kappa - 1)/(kappa + 1);
    lambda = atan(alpha1*(1-alpha2))/(1-s0); 

    z = s0 + atan(alpha1*(x-alpha2))/lambda;
    
    k = (alpha1^2*(x-alpha2).^2 + 1);
    dz = alpha1/lambda./k;
        
    %Inverse
    
    %x     = alpha2 + tan( (z - z0)*lambda )/alpha1;
    Z = lambda*(z-s0);
    dx     = lambda/alpha1 .* sec(Z).^2;

    ddx    = 2*lambda^2/alpha1 * tan(Z).*sec(Z).^2;
    
    dddx   = 2*lambda^3/alpha1 * (2*sin(Z).^2 + 1).*sec(Z).^4;
    
    ddddx  = 8*lambda^4/alpha1 * tan(Z).*sec(Z).^2.*(tan(Z).^2 + 2*sec(Z).^2);            
     
end