function [ErrAbs,ErrRel] = L2Linfty(MS1,solEx,solNum,lenx, nSpecies)

    ErrAbsVec = [];
    ErrRelVec = [];
    Int = MS1.Int;

    if lenx < length(solEx(:,1))
        disp('use vector norm')
    else
        for iSpec = 1:nSpecies
            iEx = solEx((iSpec-1)*lenx + 1:iSpec*lenx,:);
            iNum = solNum((iSpec-1)*lenx + 1:iSpec*lenx,:);
    
            ErrIAbsS = (Int*((iEx-iNum).^2)).^(1/2);
            ErrIRelS = ErrIAbsS/((Int*((iEx).^2)).^(1/2)+10^-10);
    
            ErrIAbs = max(ErrIAbsS);
            ErrIRel = max(ErrIRelS);
    
            ErrAbsVec = [ErrAbsVec,ErrIAbs];
            ErrRelVec = [ErrRelVec,ErrIRel];
        end
    end
    ErrAbs = max(abs(ErrAbsVec));
    ErrRel = max(ErrRelVec);
end