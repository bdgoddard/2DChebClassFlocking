function dataBall = Intersect_Ball(MainShape,ballShape)

    intersectFnStr = ['Intersect_Ball_' class(MainShape)];

    if(exist(intersectFnStr,'file')==2) % is a Matlab file
        intersectFn = str2func(intersectFnStr);
        area = intersectFn(ballShape,MainShape);
    else
        exc = MException('Intersect_Ball',['case ' intersectFnStr ' not implemented']);
        throw(exc);
    end
  
    y10 = ballShape.Origin(1);   
	y20 = ballShape.Origin(2);       
    
    % get polar coordinates for, e.g., FMT kernels
    dataBall.pts       = area.GetCartPts(); 
    ptsLoc.y1_kv       = dataBall.pts.y1_kv - y10;
    ptsLoc.y2_kv       = dataBall.pts.y2_kv - y20;        
    dataBall.ptsPolLoc = Cart2PolPts(ptsLoc);

    % deal with intersections at infinity
    if(y10 == inf)
        dataBall.pts.y1_kv = dataBall.pts.y1_kv + y10;
    end    
    if(y20 == inf)
        dataBall.pts.y2_kv = dataBall.pts.y2_kv + y20;
    end

    [dataBall.int,dataBall.area]     = area.ComputeIntegrationVector();                                   
    
end   