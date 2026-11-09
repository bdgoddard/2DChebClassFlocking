function res = ComputeMSMSCostFunctional(MS1,opts,rhoFull,rhoHatFull,w)

    nSpecies = opts.optsNum.nSpecies;
    lenx = length(MS1.Pts.y1_kv);
    bet = opts.optsPhys.Params.beta;

    ge.yMin = 0;
    ge.yMax = opts.optsNum.T;
    ge.N = opts.optsNum.n + 1;

    InterLine = SpectralLine(ge);   
    times = InterLine.Pts.y;

    TInt = InterLine.ComputeIntegrationVector;
    SInt = MS1.Int;
    grad = MS1.Diff.grad;

    if ~isnumeric(rhoHatFull)
        rhoHatFun = str2func(rhoHatFull);
        for iT = 1:opts.optsNum.n + 1
            t = times(iT);
            rhoHat(:,iT) = rhoHatFun(MS1,opts,t);
        end
    else
        rhoHat = rhoHatFull;
    end

    J1 = 0;

    for iS = 1:nSpecies
        rhoi = rhoFull((iS-1)*lenx+1:iS*lenx,:);
        rhoHati = rhoHat((iS-1)*lenx+1:iS*lenx,:);
        J1i = TInt*(SInt*((rhoi - rhoHati).^2))'; 
        J1 = J1 + J1i;
    end

    for iT = 1:opts.optsNum.n+1
        wt = w(:,iT);
        dProdw(:,iT)= MS1.dotProduct(wt,wt);
    end

    J2 = TInt*((SInt*(dProdw))');

    Jfun = (1/2)*J1 + (bet/2)*J2;
    
%     res = Jfun;
    res.Jfun = Jfun;
    res.J1 = J1;
    res.J2 = J2;   
end