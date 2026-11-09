function [x,Int] = FourierInt( N )
    x    =  FourierSeq(N)/2/pi;
    Int  =  ones(size(x))/N;
    Int  =  Int';
end

