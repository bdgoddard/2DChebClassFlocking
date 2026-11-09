function A = barychebevalMatrix(xIn,xOut)
% barychebevalMatrix(xIn,xOut) gives A that computes the values fOut taken
% by a rational  interpolant with barycentric weights w, collocation points
% xIn, and function values fIn at the points.
%
% Inputs: xIn and xOut, the points at which the values are known, and those
% at which it is to be interpolated. Note that xIn must be in ascending
% order.
% Outputs: A, the matrix such that fOut \approx A*fIn.

% See, e.g.,
% Barycentric Lagrange Interpolation, Jean-Paul Berrut Lloyd N. Trefethen
% SIAM REVIEW Vol. 46, No. 3, pp. 501–517

nIn = length(xIn);
nOut = length(xOut);

% weights for barycentric interpolation
% these are NOT the same as the Clenshaw-Curtis integration weights
w = (-1).^(0:(nIn-1))';
w(1) = 0.5*w(1);
w(nIn) = 0.5*w(nIn);

% sort input and output vectors, and keep track of the indexing
[sorted,indices] = sort([xIn(2:nIn);xOut(1:nOut)]);
% determine which if xIn is immediately to the left of each point
sorted(indices) = cumsum(indices<=nIn-1)+1;
index = sorted(nIn:nIn+nOut-1);

maskSame = (xIn(index)==xOut); % points in xOut where we already know the value
maskDiff = (maskSame==0);         % points in xOut ~= xIn
xDiff   = xOut(maskDiff);        % store in xDiff points s.t. xOut ~= xIn

% set up preallocated arrays
denom = zeros(length(xDiff),1);
A     = zeros(nOut,nIn);

% compute interpolant for xDiff (see paper)
if(~isempty(xDiff))
    if(length(xDiff)==1)
        temp = w./(xDiff-xIn);
        A(maskDiff,:) = temp;
        denom = denom + sum(temp);
    else
        for i=1:nIn                        % go through all points in x
            temp = w(i)./(xDiff-xIn(i));   % Paper p507: temp = c(j)./xdiff
            A(maskDiff,i) = temp;
            denom = denom + temp;
        end
    end
end

%A(maskDiff,:)  =  diag(1./denom) *  A(maskDiff,:) ;
A(maskDiff,:)  =  A(maskDiff,:) ./ denom;

% For points already in xIn we know the values
for i=1:length(maskSame)
    if(maskSame(i)==1)
        A(i,index(i)) = 1;
    end
end

end
