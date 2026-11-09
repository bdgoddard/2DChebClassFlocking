function [z,dz,dx,ddx,dddx,ddddx,ddz] = TanIntervalMap(x,alpha1,alpha2,y1,y2)

    a2  = InvLinearMap(alpha2,y1,y2);
    L = (y2-y1)/2;
    a1 = alpha1/L;
    
    [z1,dz1,dx1,ddx1,dddx1,ddddx1,ddz1] = TanMap(x,a1,a2);        
    [z,dzt]  = LinearMap(z1,y1,y2);
    dz       = dz1.*dzt;
    ddz      = ddz1.*dzt; % + ddzt .* (dz1).^2
    
    dx       = dx1/L;
    ddx      = ddx1/L.^2;
    dddx     = dddx1/L.^3;
    ddddx    = ddddx1/L.^4;
end