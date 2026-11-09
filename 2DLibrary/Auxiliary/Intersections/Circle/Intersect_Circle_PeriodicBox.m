function line = Intersect_Circle_PeriodicBox(circle,box)
    
    %Initialization
    y10      = circle.Origin(1);
    y20      = circle.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    
    R = circle.R;

    N = circle.N;
    
    central = ( (y10 >= left + R) && (y10 <= right-R) );
    
    Win   = (y10 < left + R) && (y10 >= left);
    Wout  = (y10 < left) && (y10 >= left - R);

    Ein  = (y10 <= right) && (y10 > right - R);
    Eout = (y10 <= right + R) && (y10 > right);

    edgeIn = (Ein || Win);
    edgeOut = (Eout || Wout);

    if(central)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10,y20];
        line = Circle(shape);

    elseif(edgeIn || edgeOut)
        shape.N = N;
        shape.R  = R;
        shape.Origin = [y10;y20];

        if(Ein || Eout)
            shape.h = right - y10;
            shape.WallPos = 'E';
            line = Arc(shape); 
        elseif(Win || Wout)
            shape.h = y10 - left;
            shape.WallPos = 'W';
            line = Arc(shape); 
        end
        
    else

        line = ArbitraryIntersection();
    end
end