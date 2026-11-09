function [area,Shapes] = Complement_Disc_Disc(disc1,disc2)
    %Initialization
    
    R1 = disc1.R;
    R2 = disc2.R;
    if(R1 ~= R2)
        exc = MException('Complement_Disc_Disc','discs must have same radius');
        throw(exc);
    else
        R = R1;
    end
    
    % need more points than are used in the discs?
    NFinite = disc1.N1;
    NInfinite = 2*disc1.N1;

    % N = [disc1.N1;disc1.N2];
    
    Origin1 = disc1.Origin;
    Origin2 = disc2.Origin;
    
    dx = Origin2(1) - Origin1(1);
    dy = Origin2(2) - Origin1(2);

    d = sqrt( dx^2 + dy^2 );
    
    theta = atan2(dy,dx);
    shift = Origin1;
    
    L = 5*R; % default L for infinite domains (half points are in [0,L]);

    if(d==0) % discs overlap - use an infinite annulus
        
        geomInfAnnulus.N = [NInfinite; NFinite];
        geomInfAnnulus.RMin = R;
        geomInfAnnulus.Origin = Origin1;
        geomInfAnnulus.L = L;
        Shapes(1).Shape = InfAnnulus(geomInfAnnulus);
        area = ArbitraryIntersection(Shapes);
        
    elseif(~isfinite(d)) % one or two discs at infinity
        
        isInfinite1 = any(~isfinite(Origin1)); % disc1 at infinity
        isInfinite2 = any(~isfinite(Origin2)); % disc2 at infinity
        
        if(isInfinite1 && isInfinite2) % both discs at infinity

            geomInfSpace.N = [NInfinite; NInfinite];
            geomInfSpace.L1 = L;
            geomInfSpace.L2 = L;
            Shapes(1).Shape = InfSpace(geomInfSpace);
            area = ArbitraryIntersection(Shapes);
            return; % don't rotate/shift due to infinities

        else

            geomInfAnnulus.N = [NInfinite; NFinite];
            geomInfAnnulus.RMin = R;
            geomInfAnnulus.L = L;

            if(isInfinite1) % infinite annulus around disc2
                geomInfAnnulus.Origin = Origin2;
            else            % infinite annulus around disc1
                geomInfAnnulus.Origin = Origin1;
            end

            Shapes(1).Shape = InfAnnulus(geomInfAnnulus);
            area = ArbitraryIntersection(Shapes);

            return; % don't rotate/shift
        end
    
    else % disks are both finite - need different cases if they overlap or not

        % We orient this top to bottom and then rotate because all of the
        % existing shapes have horizontal walls...
        % Discs are at 0 and d with 0 below

        % bottom half space
        geomHS_B.N = [NInfinite; NInfinite];
        geomHS_B.L1 = L;
        geomHS_B.L2 = L;
        geomHS_B.y2Min = R; % this is -R but we are going to flip the HS
        Shapes(1).Shape = HalfSpace(geomHS_B);
        Shapes(1).Shape.Pts.y2_kv = -Shapes(1).Shape.Pts.y2_kv; % flip the points in y2

        % two half-infinite strips around disc1
        geomSMD_B.N = [NInfinite; NFinite];
        geomSMD_B.R = R;
        geomSMD_B.OriginDisk = [0;0];
        geomSMD_B.L1 = L;
        geomSMD_B.y2Min = -R;
        geomSMD_B.y2Max = min(d/2,R); % if they overlap we have d/2, otherwise R

        geomSMD_BL = geomSMD_B;
        geomSMD_BL.LeftRight = 'Left';
        
        geomSMD_BR = geomSMD_B;
        geomSMD_BR.LeftRight = 'Right';

        Shapes(2).Shape = SlitMinusDisk(geomSMD_BL);
        Shapes(3).Shape = SlitMinusDisk(geomSMD_BR);

        % two half-infinite strips around disc2

        geomSMD_T.N = [NInfinite; NFinite];
        geomSMD_T.R = R;
        geomSMD_T.OriginDisk = [0;d];
        geomSMD_T.L1 = L;
        geomSMD_T.y2Min = max(d/2,d - R); % if they overlap we have d/2, otherwise d-R
        geomSMD_T.y2Max = d + R;

        geomSMD_TL = geomSMD_T;
        geomSMD_TL.LeftRight = 'Left';
        
        geomSMD_TR = geomSMD_T;
        geomSMD_TR.LeftRight = 'Right';

        Shapes(4).Shape = SlitMinusDisk(geomSMD_TL);
        Shapes(5).Shape = SlitMinusDisk(geomSMD_TR);

        % top half space
        geomHS_T.N = [NInfinite; NInfinite];
        geomHS_T.L1 = L;
        geomHS_T.L2 = L;
        geomHS_T.y2Min = d + R;
        Shapes(6).Shape = HalfSpace(geomHS_T);

        % middle infinite strip only in non-overlapping case
        if(d > 2*R)

            geomIC.N = [NInfinite; NFinite];
            geomIC.L1 = L;
            geomIC.y2Min = R;
            geomIC.y2Max = d-R;
            % Shapes(7).Shape = InfCapillary(geomIC);
            % Pick a sensible lengthscale for the finite direction - two circle radii or half the distance between their closest points (which replicates CGL points if d > 3R)
            geomIC.L2 = min(2*R, d/2 - R);
            geomICcell = namedargs2cell(geomIC);
            
            Shapes(7).Shape = InfWideCapillary(geomICcell{:});

        end

        area = ArbitraryIntersection(Shapes);
        area.RotatePts(-pi/2);
        
    end
    
    % rotate and shift to deal with spheres not centred at [0,0] and [d,0]
    area.RotatePts(theta);
    area.ShiftPts(shift)
    
end