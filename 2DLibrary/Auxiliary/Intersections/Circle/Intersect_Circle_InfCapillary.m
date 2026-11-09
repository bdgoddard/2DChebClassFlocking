function line = Intersect_Circle_InfCapillary(circle,infcapillary)
    
    %Initialization
    y10      = circle.Origin(1);
    y20      = circle.Origin(2);
    
    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end
    
    top      = infcapillary.y2Max;
    bottom   = infcapillary.y2Min;
    
    R = circle.R;
        
    N = circle.N;
         
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
        line = Circle(shape);
        
    elseif(edgeIn || edgeOut)
        shape.N = N;
        shape.R  = R;
        shape.Origin = [y10;y20];
        if(bottomIn || bottomOut)
            shape.h = y20 - bottom;
            shape.WallPos = 'S';
        elseif(topIn || topOut)
            shape.h = top - y20;
            shape.WallPos = 'N';
        end

        line = Arc(shape); 
        
    else

        line = ArbitraryIntersection();
    end
    
end