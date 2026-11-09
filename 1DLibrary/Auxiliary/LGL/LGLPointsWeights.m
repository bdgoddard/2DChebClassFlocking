function [x,w,P] = LGLPointsWeights(N)


% We have (x^2-1) P_k''(x) + 2 x P_k'(x) - k(k+1) P_k(x) = 0
% Note that [(x^2-1) P_k'(x)]' = (x^2-1) P_k''(x) + 2x P_k'(x)
% So by the above [(x^2-1) P_k'(x)]' = k(k+1) P_k(x)
% We also have Bonnet's recursion formula
% (x^2 -1) P_N'(x) = N(x P_N(x) - P_{N-1}(x))
% Then [ (x^2-1)P_N'(x) ] / [(x^2-1) P_k'(x)]' 
% = [ x P_N(x) - P_{N-1}(x) ] / [ (N+1)P_N(x) ]

% LGL are the zeros of (1-x^2) P_N'(x)

% N + 1 for shorthanding later
N1 = N+1;

% Use the Chebyshev collocation points as the first guess
% and ensure ordered from -1 to 1 with minus sign
x = -cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix, where column k is L_{k-1}(x)
P = zeros(N1,N1);

% Compute via Newton-Raphson
xold = 2;
while max(abs(x-xold))>eps
    xold = x;
        
    P(:,1) = 1;
    P(:,2 )= x;
    
    for k=2:N
        P(:,k+1) = ( (2*k-1)*x.*P(:,k) - (k-1)*P(:,k-1) )/k;
    end
     
    x = xold - ( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
             
end

% set the integration weights
w = 2./(N*N1*P(:,N1).^2);
w = w'; % which should be a row vector

end