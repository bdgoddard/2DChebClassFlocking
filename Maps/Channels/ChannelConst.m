function [y,dy,ddy] = ChannelConst(x,params)

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

y = shift*ones(size(x));
dy = zeros(size(x));
ddy = zeros(size(x));

end

