    function sol = InterpVext(MS1,BCs1,opts,t)
        InterLine = opts.optsNum.InterLine;
        VextIn = opts.optsPhys.DataInput.VextIn;
        InterpT = InterLine.ComputeInterpolationMatrix(InterLine.CompSpace(t));
        Interp=InterpT.InterPol; 
        VextOut = Interp*VextIn;

        sol = VextOut';  
    end