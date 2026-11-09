function Diff = FourierDiff(x,order)

    Diff = struct();
    if(nargin == 1)
        order = 2;
    end    

    n = length(x);      % DNS: set n as the number of discretisation points

    if(n==1)
        Diff.Dx    = 0;
        Diff.DDx   = 0;
        Diff.DDDx  = 0;
        Diff.DDDDx = 0;
        return;
    end
        
    h         = 2*pi/n;
    column    = [0 .5*(-1).^(1:n-1).*cot((1:n-1)*h/2)];
    row       = column([1 n:-1:2]);
    Diff.Dx   = toeplitz(column,row)*2*pi;% 2pi is needed to go back to the -1;1 - domain;    
    
    
    if (order == 1), return; end
 
    column2   = [(-(pi^2)/(3*h^2) - 1/6)  -.5*(-1).^(1:n-1)./((sin((1:n-1)*h/2)).^2)];
    row2      = column2([1 n:-1:2]);
    Diff.DDx  = toeplitz(column2,row2)*(2*pi)^2;
    
    if (order == 2), return; end
    
end