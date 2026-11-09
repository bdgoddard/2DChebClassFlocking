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
    top      = box.y2Max;
    bottom   = box.y2Min;
    
    L2 = abs(top-bottom);
    
    R = circle.R;

    N = circle.N;
    
    centralX = ( (y10 >= left + R) );
    centralY = ( (y20 >= bottom + R) && (y20 <= top-R) );

    topIn  = (y20 <= top) && (y20 > top - R);
    topOut = (y20 <= top + R) && (y20 > top); 

    bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R);

    leftIn   = (y10 < left + R) && (y10 >= left);
    leftOut  = (y10 < left) && (y10 >= left - R);

    containsNW = ( (y10-left)^2 + (y20-top)^2 < R^2);
    containsSW = ( (y10-left)^2 + (y20-bottom)^2 < R^2);

    Nin  = topIn && centralX;
    Nout = topOut && (centralX || (leftIn && ~containsNW) );

    Sin  = bottomIn && centralX;
    Sout = bottomOut && (centralX || (leftIn &&~containsSW) );

    Win  = leftIn && centralY;
    Wout = leftOut && (centralY || (topIn && ~containsNW) || (bottomIn && ~containsSW) );

    
    central = centralX && centralY;

    NWw = containsNW;
    SWw = containsSW;

    NWwo = leftIn && topIn && ~containsNW;
    SWwo = leftIn && bottomIn && ~containsSW;

    corner = (NWw || SWw);

    doubleCut = (NWwo || SWwo);

    edgeIn = (Nin || Sin || Win);
    edgeOut = (Nout || Sout || Wout);

    if(central)
        shape.N = N;
        shape.R = R;
        shape.Origin = [y10,y20];
        line = Circle(shape);
    elseif(corner)
        shape.N = N;
        shape.Origin = [y10;y20];
        shape.R = R;

        if(SWw)
            shape.Corner = 'SW';
            shape.CornerPos = [left;bottom];
        elseif(NWw)
            shape.Corner = 'NW';
            shape.CornerPos = [left;top];
        end

        Shapes(1).Shape = Arc(shape);
        
        % determine the periodic intersection
        shapePer = shape;
        if(SWw)
            shapePer.Corner = 'NW';
            shapePer.CornerPos = [left;top];
            shapePer.Origin = [y10;(y20-bottom)+L2];
        elseif(NWw)
            shapePer.Corner = 'SW';
            shapePer.CornerPos = [left;bottom];
            shapePer.Origin = [y10;bottom-(top-y20)];
        end
        
        Shapes(2).Shape = Arc(shapePer);
        
        line = ArbitraryIntersection(Shapes);

    elseif(doubleCut)
        shape.N = N;
        shape.Origin = [y10;y20];
        shape.R = R;

        if(SWwo)
            shape.Corner = 'SW';
            shape.CornerPos = [left;bottom];
        elseif(NWwo)
            shape.Corner = 'NW';
            shape.CornerPos = [left;top];
        end

        Shapes(1).Shape = DoubleCutCircle(shape);
        
        % determine the periodic intersection
        shapePer.N = N;
        shapePer.R = R;

        if(SWwo)
            shapePer.Origin = [y10;y20 + L2];
            shapePer.h = top - shapePer.Origin(2);
            shapePer.WallPos = 'N';
        elseif(NWwo)
            shapePer.Origin = [y10;y20 - L2];
            shapePer.h = shapePer.Origin(2) - bottom;
            shapePer.WallPos = 'S';
        end
        
        Shapes(2).Shape = Arc(shapePer);
        
        line = ArbitraryIntersection(Shapes);

    elseif(edgeIn || edgeOut)
        shape.N = N;
        shape.R  = R;
        shape.Origin = [y10;y20];
        if(Sin || Sout)
            shape.h = y20 - bottom;
            shape.WallPos = 'S';
            Shapes(1).Shape = Arc(shape);
            
            % determine the periodic intersection
            shapePer = shape;
            shapePer.h = -shape.h;
            shapePer.WallPos = 'N';
            shapePer.Origin = [y10;y20 + L2];
            Shapes(2).Shape = Arc(shapePer);
            
            line = ArbitraryIntersection(Shapes);
            
        elseif(Nin || Nout)
            shape.h = top - y20;
            shape.WallPos = 'N';
            Shapes(1).Shape = Arc(shape);
            
            % determine the periodic intersection
            shapePer = shape;
            shapePer.h = -shape.h;
            shapePer.WallPos = 'S';
            shapePer.Origin = [y10;y20 - L2];
            Shapes(2).Shape = Arc(shapePer);
            
            line = ArbitraryIntersection(Shapes);

        elseif(Win || Wout)
            shape.h = y10 - left;
            shape.WallPos = 'W';
            line = Arc(shape); 
        end
        
    else

        line = ArbitraryIntersection();
    end
end