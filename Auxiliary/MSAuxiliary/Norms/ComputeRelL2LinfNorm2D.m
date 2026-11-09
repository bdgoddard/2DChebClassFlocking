function RelNorm = ComputeRelL2LinfNorm2D(this, MNew, MOld, SInt, TInt)
MDiff = MNew - MOld;

if length(SInt) == length(MDiff) 
    SInt1 = SInt;
else
    SInt1 = [SInt,SInt];

end

L2MDiff = (SInt1*(MDiff.^2)').^(1/2);
L2MOld = (SInt1*((MOld).^2)').^(1/2);

% Take either relative or absolute error in space.
L2Rel = L2MDiff./(L2MOld+10^-10);
L2Err = min(L2Rel,L2MDiff);

RelNorm = max(L2Err); 
end