function [V,DV] = VextZero(q,opts)
    V = zeros(size(q));
    if(nargout == 2)
        DV = zeros(size(q));
    end
end