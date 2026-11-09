    function f_out = Interpf(MS1,BCs1,opts,t)
        InterLine = opts.optsNum.InterLine;
        fIn = opts.optsPhys.DataInput.fIn;
        InterpT = InterLine.ComputeInterpolationMatrix(InterLine.CompSpace(t));
        Interp=InterpT.InterPol; 
        fOut = Interp*fIn;

        f_out = fOut';
    end