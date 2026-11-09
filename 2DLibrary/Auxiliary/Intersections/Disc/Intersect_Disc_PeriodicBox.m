function area = Intersect_Disc_PeriodicBox(disc,box)
    
    %Initialization
    y10      = disc.Origin(1);
    y20      = disc.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    
    R = disc.R;
        
    N = [disc.N1;disc.N2];
    NT = disc.NT;
         
    %Compute intersection
    
    central = ( (y10 >= left + R) && (y10 <= right-R) );
    
    Win   = (y10 < left + R) && (y10 >= left);
    Wout  = (y10 < left) && (y10 >= left - R);

    Ein  = (y10 <= right) && (y10 > right - R);
    Eout = (y10 <= right + R) && (y10 > right);

    edgeIn = (Ein || Win);
    edgeOut = (Eout || Wout);
    
    if(central)
        shape.N      = N;
        shape.R      = R;
        shape.Origin = disc.Origin;
        area = Disc(shape);
        
    elseif(edgeIn)
        shape.NW = N;
        shape.NT = NT;
        shape.R  = R;
        shape.Origin = [y10;y20];                
                    
        if(Ein)
            shape.Wall_Y       = right;
            shape.Wall_VertHor = 'vertical';
            shape.Wall = 'E';
            area = BigSegment(shape); 
        elseif(Win)
            shape.Wall_Y       = left;
            shape.Wall_VertHor = 'vertical';
            shape.Wall = 'W';
            area = BigSegment(shape); 
        end
        

    elseif(edgeOut)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10;y20];
            
        if(Eout)
            shape.Wall_VertHor = 'vertical';
            shape.Wall_Y       = right;
            area = SegmentWall(shape); 
        elseif(Wout)
            shape.Wall_VertHor = 'vertical';
            shape.Wall_Y       = left; 
            area = SegmentWall(shape); 
        end
        
    else

        area = ArbitraryIntersection();
    end

end