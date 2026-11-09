function area = Intersect_Ball_PeriodicBox(ball,box)
    
    %Initialization
    y10      = ball.Origin(1);
    y20      = ball.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    
    R = ball.R;
        
    N = [ball.N1;ball.N2];
         
    %Compute intersection
    
    central = ( (y10 >= left + R) && (y10 <= right-R) );
    
    Win   = (y10 < left + R) && (y10 >= left);
    Wout  = (y10 < left) && (y10 >= left - R);

    Ein  = (y10 <= right) && (y10 > right - R);
    Eout = (y10 <= right + R) && (y10 > right);

    W = Win || Wout;
    E = Ein || Eout;
    
    shape.N      = N;
    shape.R      = R;
    shape.Origin = ball.Origin;
    
    if(central)
        shape.theta1 = 0;
        shape.theta2 = pi;   
        area  = Ball(shape); 
        
    elseif(W)
                    
        th                = acos((y10 - left)/R);
        shape.theta1      = 0;
        shape.theta2      = pi-th;        
        shape.theta0      = -pi/2;
        area  = Ball(shape);         

    elseif(E)
        th                = acos((right-y10)/R);
        shape.theta1      = th;
        shape.theta2      = pi; 
        shape.theta0      = -pi/2;
        area  = Ball(shape);        
    else

        area = ArbitraryIntersection();
    end

end