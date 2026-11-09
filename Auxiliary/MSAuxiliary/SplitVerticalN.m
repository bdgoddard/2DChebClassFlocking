function sStruct = SplitVerticalN(s,nSpecies,lenx)
    for iSpec =1:nSpecies
        siStr = strcat("s",num2str(iSpec));
        sStruct.(siStr) = s((iSpec-1)*lenx + 1: iSpec*lenx);
    end
end