function shapeParams = ShiftGeometryOrigin(shapeParams,origin)

    if(strcmp(shapeParams.shape,'Disc'))
        shapeParams.Origin = origin;
        
    elseif(strcmp(shapeParams.shape,'Box'))
        shapeParams.y1Min = shapeParams.y1Min + origin(1);
        shapeParams.y2Min = shapeParams.y2Min + origin(2);
        shapeParams.y1Max = shapeParams.y1Max + origin(1);
        shapeParams.y2Max = shapeParams.y2Max + origin(2);
        
    else
        error('ShiftGeoemetryOrigin:incorrectShape',...
       'Error. shapeParams.shape must be a valid shape');
    end

end