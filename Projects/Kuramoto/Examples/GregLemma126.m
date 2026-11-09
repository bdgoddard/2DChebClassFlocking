function GregLemma126

geom.N = 100;  % number of points to use for the pseudospectral method
geom.yMin = 0;    % boundaries of periodic domain
geom.yMax = 1;

geomPlot.N = 200;
geomPlot.yMin = geom.yMin;
geomPlot.yMax = geom.yMax;

aLine = FourierLine(geom);
aLine.ComputeAll(geomPlot);

q = aLine.Pts.y;
Int = aLine.Int;

beta = 2;

optsKernel.kernelFn = 'KuramotoKernelFourier';
optsKernel.a_n = [-1];
optsKernel.L = aLine.L;

convOpts.geom = geom;
convOpts.optsKernel = optsKernel; % note this is updated with L

saveDirBase = ['Kuramoto' filesep 'Overdamped'];
convOpts.convDir = [saveDirBase filesep 'Convolutions' filesep convOpts.optsKernel.kernelFn];
Conv = ComputeConvolutionWithParameters(aLine,convOpts);

optsVext.VextFn = 'VextCosSeries';
optsVext.a_n = -0.5;
optsVext.L = aLine.L;

VextFn = str2func(optsVext.VextFn);
V = VextFn(q,optsVext);

aList = [-4:0.1:4];
na = length(aList);
EList = zeros(size(aList));
dEList = zeros(size(aList));
for ia = 1:na
    rhoi = steadyState(aList(ia),q);
    EList(ia) = EMF(rhoi);
    dEList(ia) = max(abs(dEMF(rhoi)));
end

figure
plot(aList,EList,'r');
hold on
plot(aList,dEList,'b');

function E = EMF(rho)
    E = beta.^-1 * Int*(rho.*log(rho) + V.*rho + 1/2*rho.*(Conv*rho));
end

function dEdRho = dEMF(rho)
    dEdRho = beta.^-1 * log(rho) + 1 + V + Conv*rho;
end


function rho = steadyState(a,q)
    rho = exp(a*cos(2*pi*q));
    Z = Int*rho;
    rho = rho/Z;
end

end