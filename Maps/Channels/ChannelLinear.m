function [y,dy,ddy] = ChannelLinear(x,params)

% x = [-1,1]

shift = 0;
omega = 0;
phi   = 0;

if(isfield(params,'shift'))
    shift = params.shift;
end
if(isfield(params,'omega'))
    omega = params.omega;
end
if(isfield(params,'phi'))
    phi = params.phi;
end

y = omega*x + shift;
dy = omega*ones(size(x));
ddy = zeros(size(x));
end

