function line = Intersect_Circle_PeriodicHalfSpace(circle,box)
    % y1 spectral with wall at y1Min (left); there is no right
    % y2 periodic with top and bottom

    %Initialization
    y10      = circle.Origin(1);
    y20      = circle.Origin(2);

    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end    
        
    left     = box.y1Min;
    
    R = circle.R;

    N = circle.N;
    
    central = ( (y10 >= left + R) );

    Win   = (y10 < left + R) && (y10 >= left);
    Wout  = (y10 < left) && (y10 >= left - R);

    if(central)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10,y20];
        line = Circle(shape);

    elseif(Win || Wout)
        shape.N = N;
        shape.R  = R;
        shape.Origin = [y10;y20];
        shape.h = y10 - left;
        shape.WallPos = 'W';
        line = Arc(shape); 
        
    else

        line = ArbitraryIntersection();
    end
end