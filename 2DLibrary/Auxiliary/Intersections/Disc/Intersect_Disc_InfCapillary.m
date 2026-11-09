function area = Intersect_Disc_InfCapillary(disc,infcapillary)
    
    %Initialization
    y10      = disc.Origin(1);
    y20      = disc.Origin(2);
    
    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end
    
    top      = infcapillary.y2Max;
    bottom   = infcapillary.y2Min;
    
    R = disc.R;
        
    N = [disc.N1;disc.N2];
    NT = disc.NT;
         
    %Compute intersection
    
    central = ( (y20 >= bottom + R) && (y20 <= top-R) );
    
    topIn  = (y20 <= top) && (y20 > top - R); 
    topOut = (y20 <= top + R) && (y20 > top); 

    bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R);
    
    edgeIn = (topIn | bottomIn);
    edgeOut = (topOut | bottomOut);
    
    if(central)
        shape.N      = N;
        shape.R      = R;
        shape.Origin = [y10;y20];
        area = Disc(shape);
        
    elseif(edgeIn)
        shape.NW = N;
        shape.NT = NT;
        shape.R  = R;
        shape.Origin = [y10;y20];                
        shape.Wall_VertHor = 'horizontal';
        
        if(topIn)
            shape.Wall_Y       = top;
            shape.Wall = 'N';
        end
        if(bottomIn)
            shape.Wall_Y = bottom;
            shape.Wall = 'S';
        end
        
        area = BigSegment(shape);        

    elseif(edgeOut)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10;y20];
        shape.Wall_VertHor = 'horizontal';
        
        if(topOut)
            shape.Wall_Y       = top;         
        end
        if(bottomOut)
            shape.Wall_Y       = bottom;
        end
        
        area = SegmentWall(shape);             
        
    else

        area = ArbitraryIntersection();
    end
    
end