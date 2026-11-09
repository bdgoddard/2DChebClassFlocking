function I = GetIndicesFourier(this)

x = this.Pts.x;
y = this.Pts.y;

N = length(x);

left = false(size(x));
right = false(size(x));

bound = (left | right);

% leftFinite  = isfinite(y(left));
% rightFinite = isfinite(y(right));

%finite   = ( leftFinite & left) | ( rightFinite & right);
%infinite = (~leftFinite & left) | (~rightFinite & right);

finite = [];
infinite = [];

nLeft  = zeros(N);
nRight = zeros(N);

% nLeft(left,left) = -1;
% nRight(right,right) = 1;

% nFinite   = leftFinite*nLeft(finite,:) + rightFinite*nRight(finite,:);
% nInfinite = (~leftFinite)*nLeft(infinite,:) + (~rightFinite)*nRight(infinite,:);

nFinite = [];
nInfinite = [];

I = struct('left',left,'right',right,'bound',bound, ...
           'normalLeft',nLeft(left,:), ...
           'normalRight',nRight(right,:), ...
           'normal', nLeft(bound,:) + nRight(bound,:), ...
           'finite',finite,'infinite',infinite, ...
           'normalFinite',nFinite,'normalInfinite',nInfinite);
end