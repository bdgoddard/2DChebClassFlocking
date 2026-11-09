function dataSphere = Intersect_Sphere(MainShape,sphereShape)

    intersectFnStr = ['Intersect_Sphere_' class(MainShape)];

    if(exist(intersectFnStr,'file')==2) % is a Matlab file
        intersectFn = str2func(intersectFnStr);
        area = intersectFn(sphereShape,MainShape);
    else
        exc = MException('Intersect_Sphere',['case ' intersectFnStr ' not implemented']);
        throw(exc);
    end

    y10 = sphereShape.Origin(1);   
	y20 = sphereShape.Origin(2);       
    
    % get polar coordinates for, e.g., FMT kernels
    dataSphere.pts       = area.GetCartPts();           
    ptsLoc.y1_kv       = dataSphere.pts.y1_kv - y10;
    ptsLoc.y2_kv       = dataSphere.pts.y2_kv - y20;        
    dataSphere.ptsPolLoc = Cart2PolPts(ptsLoc);

    % deal with intersections at infinity
    if(y10 == inf)
        dataSphere.pts.y2_kv = dataSphere.pts.y2_kv + y10;
    end
    if(y20 == inf)
        dataSphere.pts.y2_kv = dataSphere.pts.y2_kv + y20;
    end

    [dataSphere.int,dataSphere.area]     = area.ComputeIntegrationVector();                                   
    
end   