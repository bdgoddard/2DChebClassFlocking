function [ErrAbs,ErrRel] = L2(MS1,solEx,solNum,lenx,nSpecies)

    ErrAbsVec = [];
    ErrRelVec = [];
    SInt = MS1.Int;
    if (lenx < length(solEx(:,1))) && (nSpecies == 1)
        disp('use vector norm')
    else
        for iSpec = 1:nSpecies
            iEx = solEx((iSpec-1)*lenx + 1:iSpec*lenx,:);
            iNum = solNum((iSpec-1)*lenx + 1:iSpec*lenx,:);
            
            ErrIAbs = (SInt*((iEx-iNum).^2)).^(1/2);
            ErrIRel = ErrIAbs/((SInt*((iEx).^2)).^(1/2) + 10^-10);
            ErrAbsVec = [ErrAbsVec,ErrIAbs];
            ErrRelVec = [ErrRelVec,ErrIRel];
        end
    end
    ErrAbs = max(ErrAbsVec);
    ErrRel = max(ErrRelVec);
end