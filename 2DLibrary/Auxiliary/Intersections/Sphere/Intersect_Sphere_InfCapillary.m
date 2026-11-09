function area = Intersect_Sphere_InfCapillary(sphere,IC)

    r   = sphere.R;
    
    if((IC.y2Max - IC.y2Min) < 2*r)
        exc = MException('Intersect_Sphere_InfCapillary','case not implemented');
        throw(exc);                
    end    
    
    N   = [sphere.N1,sphere.N2];
    y20 = sphere.Origin(2);       

    shape.N  = N;
    shape.R  = r;     
    shape.Origin = sphere.Origin;
    shape.volume = sphere.volume;

    y2Min = IC.y2Min;
    y2Max = IC.y2Max;
    if(isa(IC,'InfCapillarySkewed'))   
        y2Min = y2Min*sin(IC.alpha);
        y2Max = y2Max*sin(IC.alpha);
    end   

    %1. find points of disk fully in InfCapillary            
    if((y20 < (y2Min + r)) && (y20 >= (y2Min - r)))
        th                = acos((y20 - y2Min)/r);
        shape.theta1      = 0;
        shape.theta2      = pi-th;    
        area  = Sphere(shape);
    elseif( (y20 >= y2Min + r) && (y20 <= y2Max - r))           
        shape.theta1 = 0;
        shape.theta2 = pi;    
        area  = Sphere(shape);
    elseif((y20 > (y2Max - r)) && (y20 <= (y2Max + r)))
        th                = acos((y2Max-y20)/r);
        shape.theta1      = th;
        shape.theta2      = pi;           
        area  = Sphere(shape);
    else
        area = ArbitraryIntersection();
    end
    
    
end