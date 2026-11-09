function line = Intersect_Circle_HalfSpace(circle,halfspace)
    
    %Initialization
    y10      = circle.Origin(1);
    y20      = circle.Origin(2);
    
    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end
    if(~isfinite(y20)) 
        y20 = 0;
    end    
    
    bottom   = halfspace.y2Min;
    
    R = circle.R;
        
    N = circle.N;
         
    %Compute intersection
    
    central = ( (y20 >= bottom + R) );
    
    bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R);
    
    if(central)
        shape.N      = N;
        shape.R      = R;
        shape.Origin = [y10;y20];
        line = Circle(shape);
        
    elseif(bottomIn || bottomOut)
        shape.N = N;
        shape.R  = R;
        shape.Origin = [y10;y20];
        shape.h = y20 - bottom;
        shape.WallPos = 'S';

        line = Arc(shape); 
        
    else

        line = ArbitraryIntersection();
    end
    
end