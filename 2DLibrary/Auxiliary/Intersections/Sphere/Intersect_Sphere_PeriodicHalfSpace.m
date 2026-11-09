function area = Intersect_Sphere_PeriodicHalfSpace(sphere,halfspace)
    % y1 spectral with wall at y1Min (left); there is no right
    % y2 periodic with top and bottom

    %Initialization
    y10   = sphere.Origin(1);
    y20   = sphere.Origin(2);
    wall  = halfspace.y1Min;   
    
    % not sure if we need this
    if(~isfinite(y10)) 
        y10 = 0;
    end    
    
    left     = halfspace.y1Min;    
    
    R = sphere.R;
        
    N = [sphere.N1;sphere.N2];
         
    %Compute intersection
    
    central = ( (y10 >= left + R) );
    
    Win   = (y10 < left + R) && (y10 >= left);
    Wout  = (y10 < left) && (y10 >= left - R);

    shape.N      = N;
    shape.R      = R;
    shape.Origin = [y10;y20];
    shape.volume = sphere.volume;
    
    if(central)
        shape.theta1 = 0;
        shape.theta2 = pi;       
        area  = Sphere(shape);       
        
    elseif(Win || Wout)        
        th = acos((y10 - wall)/R);
        shape.theta1      = 0;
        shape.theta2      = pi-th;    
        shape.theta0 = -pi/2;
        area  = Sphere(shape);        
        
    else

        area = ArbitraryIntersection();
    end

end