function Diff = LGLDiff(x,P,order)
% x are the LGL points
% P is the Legendre Vandermonde Matrix, where column k is L_{k-1}(x)
% order is the order up to which to compute the differentiation matrices
% (currently unused, defaults to 2)

% See Canuto, Hussaini, Quarteroni, Zang, Spectral Methods: Fundamentals in
% Single Domains, Springer, 2006
% (2.3.28) for Dx and (2.4.36) for DDx, although we choose to just make the
% row sums zero by altering the diagonal elements

%---------------------- First order ---------------------------------------

N  = length(x) - 1;
N1 = length(x);

% compute the matrix of x_i - x_j, and pad the diagonal with 1's to avoid
% issues when i=j.
X = repmat(x,1,N1);
Xdiff = X-X'+eye(N1);

% P(:,N1) is L_N(x); replicate this column to make a square matrix
L = repmat(P(:,N1),1,N1);

% set the correct off-diagonal elements
Dx = L./(Xdiff.*L');

% zero out the diagonals
Dx(1:(N1+1):N1*N1) = 0;

% deal with two non-zero diagonal entries
Dx(1,1)   = -(N1*N)/4;
Dx(N1,N1) =  (N1*N)/4;

Diff.Dx = Dx;

%---------------------- Second order --------------------------------------

diagDx = diag(Dx);
diagDx = diagDx(:);
diagDx = repmat(diagDx,1,N1);

% set off-diagonal elements
% note the factor of 2 is missing from CHQZ (2.4.36), but is present in the
% errata of the original paper: Computers and Mathematics with Applications 38 (1999) 119
DDx = 2*Dx .* (diagDx - 1./Xdiff); 

% zero diagonal entries
DDx(1:(N1+1):N1*N1) = 0;

% set diagonal entries
DDx(1:(N1+1):N1*N1) = - sum(DDx,2)';

Diff.DDx = DDx;

end