function Diff = FDDiff(N)
    
    h=2/(N-1);
    
    D = sparse(2:N-1,[(3:N)],1/(2*h),N,N) ...
            - sparse(2:N-1,[1:(N-2)],1/(2*h),N,N);
        
	D(1,1) = -1/h; D(1,2) = 1/h;
    D(N,N) = 1/h;  D(N,N-1) = -1/h;
        
    Diff.Dx = sparse(D);
    Diff.DDx = sparse(D^2);
    
end