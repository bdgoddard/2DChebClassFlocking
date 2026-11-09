function data = Intersect(MainShape,SecondShape)

    mainClass = class(MainShape);
    
    % deal with Skewed and _FMT classes which have the same intersection
    % as their parent
    truncateList = {'Skewed','_FMT'};
    for iTruncate = 1:length(truncateList)
        if(endsWith(mainClass,truncateList{iTruncate}))
            mainClass = mainClass(1:end-length(truncateList{iTruncate}));
            break; % only need to truncate once
        end
    end
    
    secondClass = class(SecondShape);
   
    % determine intersection function
    intersectFnStr = ['Intersect_' secondClass '_' mainClass];
    noFile = false;
    
    % call it if it exists
    if(exist(intersectFnStr,'file')==2) % is a Matlab file
        intersectFn = str2func(intersectFnStr);
        area = intersectFn(SecondShape,MainShape);
    else % try with just the second shape
        intersectFnStrSecond = ['Intersect_' secondClass];
        if(exist(intersectFnStrSecond,'file')==2) % is a Matlab file
            intersectFn = str2func(intersectFnStrSecond);
            data = intersectFn(MainShape,SecondShape);
            return
        else
            noFile = true;
        end
    end
        
    if(noFile)
        disp(['\n mainClass = ' mainClass ', secondClass = ' secondClass]);
        exc = MException('Intersect','case not implemented');
        throw(exc);
    end

    y10 = SecondShape.Origin(1);   
	y20 = SecondShape.Origin(2);       

    % get polar coordinates for, e.g., FMT kernels
    data.pts       = area.GetCartPts();
    ptsLoc.y1_kv   = data.pts.y1_kv - y10;
    ptsLoc.y2_kv   = data.pts.y2_kv - y20;
    data.ptsPolLoc = Cart2PolPts(ptsLoc);

    % deal with intersections at infinity
    if(y10 == inf)
        data.pts.y1_kv = data.pts.y1_kv + y10;
    end
    if(y20 == inf)
        data.pts.y2_kv = data.pts.y2_kv + y20;
    end
    
    % compute the integration vector and length/area
    if(isa(SecondShape,'Circle'))
        [data.int,data.length] = area.ComputeIntegrationVector(); 
    else
        [data.int,data.area]   = area.ComputeIntegrationVector();  
    end
    
    
    
end