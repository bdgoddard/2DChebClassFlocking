function area = Intersect_Disc_PeriodicHalfSpace(disc,halfspace)
    % y1 spectral with wall at y1Min (left); there is no right
    % y2 periodic with top and bottom

    %Initialization
    y10      = disc.Origin(1);
    y20      = disc.Origin(2);
    
    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end    
    
    left     = halfspace.y1Min;    
    
    R = disc.R;
        
    N = [disc.N1;disc.N2];
    NT = disc.NT;
         
    %Compute intersection
    
    central = ( (y10 >= left + R) );
    
    Win   = (y10 < left + R) && (y10 >= left);
    Wout  = (y10 < left) && (y10 >= left - R);
    
    if(central)
        shape.N      = N;
        shape.R      = R;
        shape.Origin = disc.Origin;
        area = Disc(shape);        
        
    elseif(Win)
        shape.NW = N;
        shape.NT = NT;
        shape.R  = R;
        shape.Origin = [y10;y20];                
                    
        shape.Wall_Y       = left;
        shape.Wall_VertHor = 'vertical';
        shape.Wall = 'W';
        area = BigSegment(shape); 
        

    elseif(Wout)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10;y20];
        shape.Wall_VertHor = 'vertical';
        shape.Wall_Y       = left; 
        area = SegmentWall(shape); 
               
    else

        area = ArbitraryIntersection();
    end

end