function dataBox = Intersect_Box(MainShape,boxShape)
            
    if(isa(MainShape,'Box'))
        area           = Intersect_Box_Box(boxShape,MainShape);                
        
    else
        exc = MException('Intersect_Box','case not implemented');
        throw(exc);                
    end        
    
    dataBox.pts        = area.GetCartPts();    
    dataBox.int = area.ComputeIntegrationVector();
    dataBox.area = (area.y1Max - area.y1Min)*(area.y2Max - area.y2Min);
    
end   