    function g_out = Interpg(MS1,opts,t)
        InterLine = opts.optsNum.InterLine;
        gIn = opts.optsPhys.DataInput.gIn;
        InterpT = InterLine.ComputeInterpolationMatrix(InterLine.CompSpace(t));
        Interp=InterpT.InterPol;
        gOut = Interp*gIn';

        g_out = gOut';
    end