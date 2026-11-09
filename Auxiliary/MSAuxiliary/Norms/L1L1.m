function [ErrAbs,ErrRel] = L1L1(MS1,InterLine,solEx,solNum,lenx,nSpecies)

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
            
            ErrIAbsS = (SInt*abs(iEx-iNum))';
            ErrIRelS = ErrIAbsS/((SInt*abs(iEx))'+10^-10);
    
            ErrIAbs = TInt*ErrIAbsS;
            ErrIRel = TInt*ErrIRelS;
    
            ErrAbsVec = [ErrAbsVec,ErrIAbs];
            ErrRelVec = [ErrRelVec,ErrIRel];
        end
    end
    ErrAbs = max(ErrAbsVec);
    ErrRel = max(ErrRelVec);
end