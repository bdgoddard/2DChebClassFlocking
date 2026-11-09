function [ErrAbs,ErrRel] = L1(MS1,solEx,solNum,lenx,nSpecies)

    ErrAbsVec = [];
    ErrRelVec = [];
    SInt = MS1.Int;
    if lenx < length(solEx(:,1))
        disp('use vector norm')
    else
        for iSpec = 1:nSpecies
            iEx = solEx((iSpec-1)*lenx + 1:iSpec*lenx,:);
            iNum = solNum((iSpec-1)*lenx + 1:iSpec*lenx,:);
            
            ErrIAbs = SInt*(iEx-iNum);
            ErrIRel = ErrIAbs/(SInt*(iEx) + 10^-10);
            ErrAbsVec = [ErrAbsVec,ErrIAbs];
            ErrRelVec = [ErrRelVec,ErrIRel];
        end
    end
    ErrAbs = max(ErrAbsVec);
    ErrRel = max(ErrRelVec);
end