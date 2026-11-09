function area = Intersect_Ball_HalfSpace(ball,HS)

    r   = ball.R;
    N   = [ball.N1,ball.N2];
	y20 = ball.Origin(2);       
    
	shape.N  = N;
    shape.R  = r;     
    shape.Origin = ball.Origin;

    y2Min = HS.y2Min;
    if(isa(HS,'HalfSpaceSkewed'))   
        y2Min = y2Min*sin(HS.alpha);
    end

    %1. find points of disk in HalfSpace            
    if(y20 == inf)
        shape.Origin(2) = 0;
        shape.theta1 = 0;
        shape.theta2 = pi;
        area  = Ball(shape); 
    elseif(y20 >= y2Min + r)
        %1a. if full disk is in HalfSpace                
        shape.theta1 = 0;
        shape.theta2 = pi;   
        area  = Ball(shape); 
    elseif((y20 < (y2Min + r)) && ...
           (y20 >= (y2Min - r)))
        %1b. if part of disk is in HalfSpace  (>= half)
        %1b1. Integrate over segment in HalfSpace
        %shape.Origin    = [0,y20];
        th                = acos((y20 - y2Min)/r);
        shape.theta1      = 0;
        shape.theta2      = pi-th;       
        area  = Ball(shape); 
    else
        area = ArbitraryIntersection();
    end
    
    
    
end