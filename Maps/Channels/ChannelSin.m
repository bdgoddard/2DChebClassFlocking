function [y,dy,ddy] = ChannelSin(x,params)

% x = [-1,1]

shift = 0;
omega = 0;
phi   = 0;
A     = 1;

if(isfield(params,'shift'))
    shift = params.shift;
end
if(isfield(params,'omega'))
    omega = params.omega;
end
if(isfield(params,'phi'))
    phi = params.phi;
end
if(isfield(params,'A'))
    A = params.A;
end


y = A*sin(omega*x + phi) + shift;
dy = A*omega*cos(omega*x + phi);

ddy = -A*omega^2*sin(omega*x + phi);


end

