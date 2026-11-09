function area = Intersect_Ball_InfCapillary(ball,IC)

    r   = ball.R;
    
    if((IC.y2Max - IC.y2Min) < 2*r)
        exc = MException('Intersect_Ball_InfCapillary','case not implemented');
        throw(exc);                
    end
    
    N   = [ball.N1,ball.N2];
	y20 = ball.Origin(2);       
    
	shape.N  = N;
    shape.R  = r;     
    shape.Origin = ball.Origin;

    y2Min = IC.y2Min;
    y2Max = IC.y2Max;
    if(isa(IC,'InfCapillarySkewed'))   
        y2Min = y2Min*sin(IC.alpha);
        y2Max = y2Max*sin(IC.alpha);
    end   

    %1. find points of disk fully in InfCapillary            
    if((y20 < (y2Min + r)) && (y20 >= (y2Min - r))) %bottom
        th                = acos((y20 - y2Min)/r);
        shape.theta1      = 0;
        shape.theta2      = pi-th;        
        area  = Ball(shape); 
    elseif( (y20 >= y2Min + r) && (y20 <= y2Max - r))           
        shape.theta1 = 0;
        shape.theta2 = pi;   
        area  = Ball(shape); 
    elseif((y20 > (y2Max - r)) && (y20 <= (y2Max + r))) %top
        th                = acos((y2Max-y20)/r);
        shape.theta1      = th;
        shape.theta2      = pi;           
        area  = Ball(shape); 
    else
        area = ArbitraryIntersection();
    end
    
    
end