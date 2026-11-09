function default = KuramotoDefaultsPhaseSpace()

defaultOverdamped = KuramotoDefaults();

geom.N = [40;30];
pMax = 7;
geom.y1Min = -pMax;
geom.y1Max = pMax;
geom.y2Min = 0;
geom.y2Max = 2*pi;

pMaxPlot = 4;
geomPlot = geom;
geomPlot.y1Min = -pMaxPlot;
geomPlot.y1Max = pMaxPlot;
geomPlot.N = [50;50];

default.geom = geom;
default.geomPlot = geomPlot;

default.gamma = 1;
default.beta = defaultOverdamped.beta;
default.kappa = defaultOverdamped.kappa;

default.relTol = 10^-8;
default.absTol = 10^-8;

default.toEq = true;
default.eqTol = 10^-7;

%default.tMax = 1000;
default.tMax = defaultOverdamped.tMax;
default.nSaves = defaultOverdamped.nSaves;

default.optsIC = defaultOverdamped.optsIC;
default.optsIC.phiICFn = 'KuramotoMB';
default.optsIC.beta = default.beta;

default.optsKernel = defaultOverdamped.optsKernel;

default.optsVext = defaultOverdamped.optsVext;

default.optsSDE = defaultOverdamped.optsSDE;

end