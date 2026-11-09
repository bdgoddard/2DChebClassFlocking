function line = Intersect_Circle_PeriodicBox(circle,box)
    
    %Initialization
    y10      = circle.Origin(1);
    y20      = circle.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    top      = box.y2Max;
    bottom   = box.y2Min;
    
    L2 = abs(top-bottom);
    
    R = circle.R;

    N = circle.N;
    
    centralX = ( (y10 >= left + R) && (y10 <= right-R) );
    centralY = ( (y20 >= bottom + R) && (y20 <= top-R) );

    topIn  = (y20 <= top) && (y20 > top - R);
    topOut = (y20 <= top + R) && (y20 > top); 

    bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R);

    leftIn   = (y10 < left + R) && (y10 >= left);
    leftOut  = (y10 < left) && (y10 >= left - R);

    rightIn  = (y10 <= right) && (y10 > right - R);
    rightOut = (y10 <= right + R) && (y10 > right);

    containsNW = ( (y10-left)^2 + (y20-top)^2 < R^2);
    containsNE = ( (y10-right)^2 + (y20-top)^2 < R^2);
    containsSE = ( (y10-right)^2 + (y20-bottom)^2 < R^2);
    containsSW = ( (y10-left)^2 + (y20-bottom)^2 < R^2);


    Nin  = topIn && centralX;
    Nout = topOut && (centralX || (leftIn && ~containsNW) || (rightIn && ~containsNE) );

    Sin  = bottomIn && centralX;
    Sout = bottomOut && (centralX || (leftIn &&~containsSW) || (rightIn && ~containsSE) );

    Ein  = rightIn && centralY;
    Eout = rightOut && (centralY || (topIn && ~containsNE) || (bottomIn && ~containsSE) );

    Win  = leftIn && centralY;
    Wout = leftOut && (centralY || (topIn && ~containsNW) || (bottomIn && ~containsSW) );

    
    central = centralX && centralY;

    NWw = containsNW;
    NEw = containsNE;
    SWw = containsSW;
    SEw = containsSE;

    NWwo = leftIn && topIn && ~containsNW;
    NEwo = rightIn && topIn && ~containsNE;
    SWwo = leftIn && bottomIn && ~containsSW;
    SEwo = rightIn && bottomIn && ~containsSE;

    corner = (NWw || NEw || SEw || SWw);

    doubleCut = (NWwo || NEwo || SEwo || SWwo);

    edgeIn = (Nin || Sin || Ein || Win);
    edgeOut = (Nout || Sout || Eout || Wout);

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
        elseif(NEw)
            shape.Corner = 'NE';
            shape.CornerPos = [right;top];
        elseif(SEw)
            shape.Corner = 'SE';
            shape.CornerPos = [right;bottom];
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
        elseif(NEw)
            shapePer.Corner = 'SE';
            shapePer.CornerPos = [right;bottom];
            shapePer.Origin = [y10;bottom-(top-y20)];
        elseif(SEw)
            shapePer.Corner = 'NE';
            shapePer.CornerPos = [right;top];
            shapePer.Origin = [y10;(y20-bottom)+L2];
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
        elseif(NEwo)
            shape.Corner = 'NE';
            shape.CornerPos = [right;top];
        elseif(SEwo)
            shape.Corner = 'SE';
            shape.CornerPos = [right;bottom];
        end

        Shapes(1).Shape = DoubleCutCircle(shape);
        
        % determine the periodic intersection
        shapePer.N = N;
        shapePer.R = R;

        if(SWwo || SEwo)
            shapePer.Origin = [y10;y20 + L2];
            shapePer.h = top - shapePer.Origin(2);
            shapePer.WallPos = 'N';
        elseif(NWwo || NEwo)
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

        elseif(Ein || Eout)
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