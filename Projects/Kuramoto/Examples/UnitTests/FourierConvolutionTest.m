function FourierConvolutionTest()

geom.N = 100;  % number of points to use for the pseudospectral method
geom.yMin = 0;    % boundaries of periodic domain
geom.yMax = 1;

geomPlot.N = 200;
geomPlot.yMin = geom.yMin;
geomPlot.yMax = geom.yMax;

aLine = FourierLine(geom);
aLine.ComputeAll(geomPlot);

optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = [-1];
optsKernel.L = aLine.L;

Conv = aLine.ComputeConvolutionMatrix(@kernelWrapper);

% check of an exact formula with a_n = -1 (-cos(q))
q = aLine.Pts.y;
f = testFn(q);
optsPlot.plain = true;
optsPlot.linewidth = 1.5;
optsPlot.linecolor = 'r';
aLine.plot(Conv*f,optsPlot)
hold on
optsPlot.linecolor = 'b';
optsPlot.linestyle = '--';
aLine.plot(convExact(q),optsPlot)
max(abs(Conv*f-convExact(q)))


    function K = kernelWrapper(y)
        kernelFn = str2func(optsKernel.kernelFn);
        if(nargin(kernelFn) == 1)
            K = kernelFn(y);
        else
            K = kernelFn(y,optsKernel);
        end
    end

    function c = convExact(q)
        alpha = 2*pi/optsKernel.L;
        c = - 3/4*pi/alpha*sin(alpha*q);
    end

    function f = testFn(q)
        alpha = 2*pi/optsKernel.L;
        f = sin(alpha*q).^3;
    end
end