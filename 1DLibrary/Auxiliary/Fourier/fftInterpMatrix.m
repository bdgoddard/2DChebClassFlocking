function  Inter  = fftInterpMatrix( x, interp )
%%FFT Interpolation
LL        = 2*pi;
phi       = interp*LL;
N         = length(x);
M     = floor((N+1)/2);

n1   = 1:(M-1);       % first half of points, excluding zero
n2   = (M+1):(N-1);   % second half

interpol  = [exp(2*pi*1i*phi/LL*0) ,...
    exp(2*pi*1i*phi/LL*n1) , ...
    cos(2*pi*M*phi/LL) ,...
    exp(2*pi*1i*phi/LL*(n2-N))]/N;

n = (0:N-1);
k = (0:N-1);

FFTMatrix = exp( -1i * 2*pi *n'*k /N);
Inter = real(interpol*FFTMatrix);
end

