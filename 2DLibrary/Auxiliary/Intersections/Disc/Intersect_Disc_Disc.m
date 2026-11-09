function [area,complement,icase] = Intersect_Disc_Disc(disc1,disc2)
% N.B. The complement is that portion of *disc2* with area i.e. intersection removed

% Output argument icase defines which of the possible intersection cases occurred:
% icase = 0: discs coincide (d=0)
% icase = 1: discs intersect, and the centre of one disc is inside the other (0<d<R)
% icase = 2: discs intersect, and the centre of one disc is outside the other (R<=d<2R)
% icase = 3: discs do not intersect (2R <= d <= inf)
% icase = 4: one or both discs at infinity (d = inf)

%Initialization

R1 = disc1.R;
R2 = disc2.R;
if(R1 ~= R2)
    exc = MException('Intersect_Disc_Disc','discs must have same radius');
    throw(exc);
else
    R = R1;
end

N = [disc1.N1;disc1.N2];

Origin1 = disc1.Origin;
Origin2 = disc2.Origin;

dx = Origin2(1) - Origin1(1);
dy = Origin2(2) - Origin1(2);

d = sqrt( dx^2 + dy^2 );

theta = atan2(dy,dx);
shift = Origin1;

if(d==0)
    
    geomDisc.R = R;
    geomDisc.N = N;
    geomDisc.Origin = [0;0];
    Shapes(1).Shape = Disc(geomDisc);
    area = ArbitraryIntersection(Shapes);
    
%elseif(d==inf)
elseif(~isfinite(d)) % deal with inf and nan
    
    area = ArbitraryIntersection();
    
elseif(d<2*R)
    
    % left segment
    geomSegmentL.N = N;
    geomSegmentL.R = R;
    geomSegmentL.Wall_VertHor = 'vertical';
    geomSegmentL.Wall_Y = d/2;
    geomSegmentL.Origin = [d;0];
    Shapes(1).Shape = SegmentWall(geomSegmentL);
    
    % right segment
    geomSegmentR = geomSegmentL;
    geomSegmentR.Origin = [0;0];
    Shapes(2).Shape = SegmentWall(geomSegmentR);
    
    area = ArbitraryIntersection(Shapes);
    
else
    
    area = ArbitraryIntersection();
end

% rotate and shift to deal with spheres not centred at [0,0] and [d,0]
area.RotatePts(theta);
area.ShiftPts(shift);

% clear Shapes from earlier
Shapes = struct;

if(d==0)
    
    complement = ArbitraryIntersection();
    
%elseif(d==inf)
elseif(~isfinite(d))
    
    complement = ArbitraryIntersection();
    
elseif(d<R)
    
    % moon
    shape.N = N;
    shape.R = R;
    shape.d = d;
    Shapes(1).Shape = Moon(shape);
    complement = ArbitraryIntersection(Shapes);
    
elseif(d<2*R)
    
    % moon slice
    geomMoonSlice.N = N;
    geomMoonSlice.R = R;
    geomMoonSlice.d = d;
    Shapes(1).Shape = MoonSlice(geomMoonSlice);
    
    % top segment
    geomSegmentT.N = N;
    geomSegmentT.R = R;
    geomSegmentT.Wall_VertHor = 'horizontal';
    geomSegmentT.Wall_Y = sqrt(R^2-(d/2)^2);
    geomSegmentT.Origin = [d;0];
    Shapes(2).Shape = SegmentWall(geomSegmentT);
    
    % bottom segment
    geomSegmentB = geomSegmentT;
    geomSegmentB.Wall_Y = -geomSegmentT.Wall_Y;
    
    Shapes(3).Shape = SegmentWall(geomSegmentB);
    
    complement = ArbitraryIntersection(Shapes);
else
    % don't use disc2 directly and modify the origin as this also
    % modifies the origin of the original disc2!
    geomDisc.R = R;
    geomDisc.N = N;
    geomDisc.Origin = [d;0];
    Shapes(1).Shape = Disc(geomDisc);
    complement = ArbitraryIntersection(Shapes);
    
end

% rotate and shift to deal with discs not centred at [0,0] and [d,0]
complement.RotatePts(theta);
complement.ShiftPts(shift); % shift to origin of first disc

% Logic to define icase output argument
if(d==0)
    icase = 0;
elseif(d<R)
    icase = 1;
elseif(d<2*R)
    icase = 2;
elseif(d==inf)
    icase = 4;
else
    icase = 3;
end
icase = int32(icase);
end