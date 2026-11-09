function area = Intersect_Disc_PeriodicBox(disc,box)
    
    %Initialization
    y10      = disc.Origin(1);
    y20      = disc.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    top      = box.y2Max;
    bottom   = box.y2Min;

    L2 = abs(top-bottom);    
    
    R = disc.R;
        
    N = [disc.N1;disc.N2];
    NT = disc.NT;
         
    %Compute intersection
    
    centralX = ( (y10 >= left + R) && (y10 <= right-R) );
    centralY = ( (y20 >= bottom + R) && (y20 <= top-R) );

    topOn = (abs(y20 - top) < eps);
    topIn  = (y20 <= top) && (y20 > top - R) && ~topOn; 
    topOut = (y20 <= top + R) && (y20 > top) && ~topOn; 

    bottomOn = (abs(y20 - bottom) < eps);
    bottomIn  = (y20 < bottom + R) && (y20 >= bottom) && ~bottomOn;
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R) && ~bottomOn;
    
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
    edgeOn = (topOn || bottomOn);
    
    if(central)
        shape.N      = N;
        shape.R      = R;
        shape.Origin = disc.Origin;
        area = Disc(shape);

    elseif(corner)
        shape.NT = NT;
        shape.NS = N;
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
        
        Shapes(1).Shape = CornerDisc(shape);         
        
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
        
        Shapes(2).Shape = CornerDisc(shapePer);                
        area = ArbitraryIntersection(Shapes);
        
    elseif(doubleCut)
        shape.NT = NT;
        shape.NW = N;

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

        Shapes(1).Shape = DoubleCutDisc(shape);


        % determine the periodic intersection
        shapePer.N = N;
        shapePer.R = R;

        if(SWwo)
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall_Y       = top;
            shapePer.Origin = [y10;(y20-bottom)+L2];
        elseif(NWwo)
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall_Y       = bottom;
            shapePer.Origin = [y10;bottom-(top-y20)];
        elseif(NEwo)
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall_Y       = bottom;
            shapePer.Origin = [y10;bottom-(top-y20)];
        elseif(SEwo)
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall_Y       = top;
            shapePer.Origin = [y10;(y20-bottom)+L2];
        end
        
        Shapes(2).Shape = SegmentWall(shapePer);
        
        area = ArbitraryIntersection(Shapes);
        
        
    elseif(edgeIn)
        shape.NW = N;
        shape.NT = NT;
        shape.R  = R;
        shape.Origin = [y10;y20];                
        
        if(Sin)
            shape.Wall_Y = bottom;
            shape.Wall_VertHor = 'horizontal';
            shape.Wall = 'S';
            Shapes(1).Shape = BigSegment(shape);
            
            % periodic intersection
            shapePer.N = shape.NW;
            shapePer.R = shape.R;
            shapePer.Origin = [shape.Origin(1),shape.Origin(2)+L2];
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall_Y       = top;
            Shapes(2).Shape = SegmentWall(shapePer);
            
            area = ArbitraryIntersection(Shapes);
                        
        elseif(Nin)
            shape.Wall_Y       = top;
            shape.Wall_VertHor = 'horizontal';
            shape.Wall = 'N';
            
            Shapes(1).Shape = BigSegment(shape);
            
            % periodic intersection
            shapePer.N = shape.NW;
            shapePer.R = shape.R;
            shapePer.Origin = [shape.Origin(1),shape.Origin(2)-L2];
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall_Y       = bottom;
            Shapes(2).Shape = SegmentWall(shapePer);
            
            area = ArbitraryIntersection(Shapes);
            
        elseif(Ein)
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
        if(Sout)
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = bottom;     
            Shapes(1).Shape = SegmentWall(shape); 
            
            % periodic intersection
            shapePer.NW = N;
            shapePer.NT = NT;
            shapePer.R = R;
            shapePer.Origin = [shape.Origin(1),shape.Origin(2)+L2];
            shapePer.Wall_Y = top;
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall = 'N';
            Shapes(2).Shape = BigSegment(shapePer);

            area = ArbitraryIntersection(Shapes);
            
        elseif(Nout)
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = top;         
            Shapes(1).Shape = SegmentWall(shape); 
            
            % periodic intersection
            shapePer.NW = N;
            shapePer.NT = NT;
            shapePer.R = R;
            shapePer.Origin = [shape.Origin(1),shape.Origin(2)-L2];
            shapePer.Wall_Y = bottom;
            shapePer.Wall_VertHor = 'horizontal';
            shapePer.Wall = 'S';
            Shapes(2).Shape = BigSegment(shapePer);

            area = ArbitraryIntersection(Shapes);
            
        elseif(Eout)
            shape.Wall_VertHor = 'vertical';
            shape.Wall_Y       = right;
            area = SegmentWall(shape); 
        elseif(Wout)
            shape.Wall_VertHor = 'vertical';
            shape.Wall_Y       = left; 
            area = SegmentWall(shape); 
        end
        
    elseif(edgeOn)
        shape.NW = N;
        shape.NT = NT;
        shape.R  = R;            
        
        if(bottomOn)
            shape.Wall_Y = bottom;
            shape.Origin = [y10;bottom+eps];
            shape.Wall_VertHor = 'horizontal';
            shape.Wall = 'S';
            Shapes(1).Shape = BigSegment(shape);
            
            % periodic intersection
            shapePer = shape;
            shapePer.Wall_Y    = top;
            shapePer.Wall = 'N';
            shapePer.Origin = [shape.Origin(1),top-eps];
            Shapes(2).Shape = BigSegment(shapePer);            
            
            area = ArbitraryIntersection(Shapes);
            
        elseif(topOn)
            shape.Wall_Y = top;
            shape.Origin = [y10;top-eps];
            shape.Wall_VertHor = 'horizontal';
            shape.Wall = 'N';
            Shapes(1).Shape = BigSegment(shape);
                    
            % periodic intersection
            shapePer = shape;
            shapePer.Wall_Y    = bottom;
            shapePer.Wall = 'S';
            shapePer.Origin = [shape.Origin(1),bottom+eps];
            Shapes(2).Shape = BigSegment(shapePer);

            area = ArbitraryIntersection(Shapes);
        end
        
        
    else

        area = ArbitraryIntersection();
    end

end