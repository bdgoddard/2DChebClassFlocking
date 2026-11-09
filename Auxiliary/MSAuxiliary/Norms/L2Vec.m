function [ErrAbs,ErrRel] = L2Vec(MS1,solEx,solNum,lenx,nSpecies)

    ErrAbsVec = [];
    ErrRelVec = [];
    SInt = MS1.Int;

    if nSpecies > 1
        disp('implement!!')
    else
        Diff = solNum- solEx;
        ErrIAbs = (SInt*MS1.dotProduct(Diff,Diff)).^(1/2);
        ErrIRel = ErrIAbs/((SInt*MS1.dotProduct(solEx,solEx)).^(1/2) + 10^-10);
        ErrAbsVec = [ErrAbsVec,ErrIAbs];
        ErrRelVec = [ErrRelVec,ErrIRel];
    end
    ErrAbs = max(ErrAbsVec);
    ErrRel = max(ErrRelVec);
end