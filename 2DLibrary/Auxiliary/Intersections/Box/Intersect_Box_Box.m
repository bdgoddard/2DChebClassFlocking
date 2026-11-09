function area = Intersect_Box_Box(subBox,mainBox)
    
%     origin = subBox.Origin;
%     y10 = origin(1);
%     y20 = origin(2);

    %Initialization
    mainLeft     = mainBox.y1Min;
    mainRight    = mainBox.y1Max;
    mainTop      = mainBox.y2Max;
    mainBottom   = mainBox.y2Min;
    
    subLeft     = subBox.y1Min;
    subRight    = subBox.y1Max;
    subTop      = subBox.y2Max;
    subBottom   = subBox.y2Min;
    
    % check there is some overlap
    overlap =    (subLeft < mainRight) ...
              && (subRight > mainLeft) ...
              && (subBottom < mainTop) ...
              && (subTop > mainBottom);
          
	% create intersected box
    if (overlap)
        shape.y1Min = max(subLeft,mainLeft);
        shape.y1Max = min(subRight,mainRight);
        shape.y2Min = max(subBottom,mainBottom);
        shape.y2Max = min(subTop,mainTop);

        shape.N = [subBox.N1;subBox.N2];
        %shape.Origin = subBox.Origin;
        
        area = Box(shape);
        
    else
        area = [];
    end

end