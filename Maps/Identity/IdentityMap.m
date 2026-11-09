function  [z,dz,dx,ddx,dddx,ddddx] = IdentityMap( xf,a,b )
        z = xf;        
        dz   = ones(size(xf));
        if (nargout == 2), return; end
        
        dx   = dz;
        ddx  = zeros(size(z));
        dddx = zeros(size(z));
        ddddx = zeros(size(z));
end

