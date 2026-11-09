function [ErrAbs,ErrRel] = L2L1(MS1,InterLine,solEx,solNum,lenx,nSpecies)

    ErrAbsVec = [];
    ErrRelVec = [];
    SInt = MS1.Int;
    TInt = InterLine.ComputeIntegrationVector;
    if lenx < length(solEx(:,1))
        disp('use vector norm')
    else
        for iSpec = 1:nSpecies
            iEx = solEx((iSpec-1)*lenx + 1:iSpec*lenx,:);
            iNum = solNum((iSpec-1)*lenx + 1:iSpec*lenx,:);
         
            ErrIAbsS = ((SInt*((iEx-iNum).^2)).^(1/2))';
            ErrIRelS = ErrIAbsS/((((SInt*((iEx).^2)).^(1/2))')+10^-10);
            
            ErrIAbs = TInt*ErrIAbsS;
            ErrIRel = TInt*ErrIRelS;
    
            ErrAbsVec = [ErrAbsVec,ErrIAbs];
            ErrRelVec = [ErrRelVec,ErrIRel];
        end
    end
    ErrAbs = max(ErrAbsVec);
    ErrRel = max(ErrRelVec);
end