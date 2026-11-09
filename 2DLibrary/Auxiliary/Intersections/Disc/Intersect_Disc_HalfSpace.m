function area = Intersect_Disc_HalfSpace(disc,halfspace)
    
    %Initialization
    y10      = disc.Origin(1);
    y20      = disc.Origin(2);
    
    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end
    if(~isfinite(y20)) 
        y20 = 0;
    end

    bottom   = halfspace.y2Min;
    
    R = disc.R;
        
    N = [disc.N1;disc.N2];
    NT = disc.NT;
         
    %Compute intersection
    
    central = (y20 >= bottom + R);
    
    bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R);
        
    if(central)
        shape.N      = N;
        shape.R      = R;
        shape.Origin = [y10;y20];
        area = Disc(shape);
        
    elseif(bottomIn)
        shape.NW = N;
        shape.NT = NT;
        shape.R  = R;
        shape.Origin = [y10;y20];                
        shape.Wall_VertHor = 'horizontal';
        shape.Wall_Y = bottom;
        shape.Wall = 'S';
        
        area = BigSegment(shape);        

    elseif(bottomOut)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10;y20];
        shape.Wall_VertHor = 'horizontal';
        shape.Wall_Y       = bottom;
        
        area = SegmentWall(shape);             
        
    else

        area = ArbitraryIntersection();
    end

end