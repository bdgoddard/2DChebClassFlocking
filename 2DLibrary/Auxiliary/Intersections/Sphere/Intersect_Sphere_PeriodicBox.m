function area = Intersect_Sphere_PeriodicBox(sphere,box)
    
    %Initialization
    y10      = sphere.Origin(1);
    y20      = sphere.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    
    R = sphere.R;
        
    N = [sphere.N1;sphere.N2];
         
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
    shape.Origin = sphere.Origin;
    shape.volume = sphere.volume;
    
    if(central)
        shape.theta1 = 0;
        shape.theta2 = pi;   
        area  = Sphere(shape); 
        
    elseif(W)
                    
        th                = acos((y10 - left)/R);
        shape.theta1      = 0;
        shape.theta2      = pi-th;        
        shape.theta0      = -pi/2;
        area  = Sphere(shape);         

    elseif(E)
        th                = acos((right-y10)/R);
        shape.theta1      = th;
        shape.theta2      = pi; 
        shape.theta0      = -pi/2;
        area  = Sphere(shape);        
    else

        area = ArbitraryIntersection();
    end

end