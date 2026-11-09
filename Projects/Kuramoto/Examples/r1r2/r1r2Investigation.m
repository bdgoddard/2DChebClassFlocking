function r1r2Investigation()

geom.yMin = 0;
geom.yMax = 2*pi;
geom.N = 50;
aLine = FourierLine(geom);
aLine.ComputeAll();

Int = aLine.Int;
x = aLine.Pts.y;

K = 2.2/2;

geomR.yMin = -1;
geomR.yMax = 1;
geomR.N = 20;

geomPlotR = geomR;
geomPlotR.N = 100;

aLineR = SpectralLine(geomR);
aLineR.ComputeAll(geomPlotR);


r1 = aLineR.Pts.y;
integrals = zeros(size(r1));
integralscos = zeros(size(r1));

% K = 2.2/2
r1Sol = -0.6580;
r2Sol = 0.3712;

F = Int*(integrand(x,r1Sol,r2Sol));
Fcos = Int*(integrandcos(x,r1Sol,r2Sol));
Fcos2 = Int*(integrandcos2(x,r1Sol,r2Sol));

% Fcos/F
% Fcos2/F

r2 = r2Sol;

for iR1 = 1:length(r1)
    integrals(iR1) = Int*(integrand(x,r1(iR1),r2));
    integralscos(iR1) = Int*(integrandcos(x,r1(iR1),r2));
end

Dr1 = aLineR.Diff.Dy;

Dintegrals = Dr1*integrals;
DDintegrals = Dr1*Dintegrals;
% aLineR.plot(2*K*r1.*integrals - Dintegrals);

InterR1 = aLineR.ComputeInterpolationMatrixPhys(r1Sol);
InterPolR1 = InterR1.InterPol;

id1 =2*K*r1.*integrals - Dintegrals;
InterPolR1*id1

id2 = (2*K)^2/2*(r2+1).*integrals - DDintegrals;
InterPolR1*id2

% %id3 = DDintegrals - (2*K*integrals + (2*K*r1).^2.*integrals);
% id3 = DDintegrals - Dr1*(2*K*r1.*integrals);
% InterPolR1*id3
% 
% aLineR.plot(DDintegrals);
% hold on
% aLineR.plot( Dr1*(2*K*r1.*integrals));

%aLineR.plot(r1List.*integrals - integralscos);

%aLineR.plot(Dr1*integrals - 2*K*r1.*integrals);


    function f = integrand(x,r1,r2)
        f = exp(2*K*(r1*cos(x) + (r2/2)*cos(2*x)));
    end

    function f = integrandcos(x,r1,r2)
        f = cos(x).*integrand(x,r1,r2);
    end

    function f = integrandcos2(x,r1,r2)
        f = cos(2*x).*integrand(x,r1,r2);
    end


end