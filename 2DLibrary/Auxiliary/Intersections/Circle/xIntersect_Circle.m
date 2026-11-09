function [dataCircle] = Intersect_Circle(MainShape,circle)            

    intersectFnStr = ['Intersect_Circle_' class(MainShape)];

    if(exist(intersectFnStr,'file')==2) % is a Matlab file
        intersectFn = str2func(intersectFnStr);
        line = intersectFn(circle,MainShape);
    else
        exc = MException('IntersectCircle',['case ' intersectFnStr ' not implemented']);
        throw(exc);
    end
      
    y10     = circle.Origin(1);
    y20     = circle.Origin(2);

    % get polar coordinates for, e.g., FMT kernels
    dataCircle.pts = line.GetCartPts();
    ptsLoc.y1_kv = dataCircle.pts.y1_kv - y10;
    ptsLoc.y2_kv = dataCircle.pts.y2_kv - y20;
    dataCircle.ptsPolLoc = Cart2PolPts(ptsLoc);
        
    % deal with intersections at infinity
    if(y10 == inf)
        dataCircle.pts.y1_kv = dataCircle.pts.y1_kv + y10;
    end
    if(y20 == inf)
        dataCircle.pts.y2_kv = dataCircle.pts.y2_kv + y20;
    end

    [dataCircle.int,dataCircle.length] = line.ComputeIntegrationVector(); 
    
end

%     if(isa(MainShape,'HalfSpace'))
% 
%         % not sure what this does for skewed half spaces?
%         line = Intersect_Circle_HalfSpace(circle,MainShape);
% 
%     elseif(isa(MainShape,'Box'))
% 
%         line = Intersect_Circle_Box(circle,MainShape);
% 
%     elseif(isa(MainShape,'PeriodicBox'))
% 
%         line = Intersect_Circle_PeriodicBox(circle,MainShape);        
% 
%     elseif(isa(MainShape,'InfCapillary') && (MainShape.y2Max - MainShape.y2Min)>= 2*circle.R)   
% 
%         line = Intersect_Circle_InfCapillary(circle,MainShape);        
% 
%     else
%         exc = MException('HalfSpace_FMT:Intersect_Circle','case not implemented');
%         throw(exc);                
%     end
