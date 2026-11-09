function sMat = VerticaltoMatrix(s,nSpecies,lenx)
    sMat = zeros(lenx,nSpecies);
    for iSpec =1:nSpecies
        sMat(:,iSpec) = s((iSpec-1)*lenx + 1: iSpec*lenx);
    end
end