function Interpol = FDInterpolation(x,xPlot)

nPlot = length(xPlot);
N = length(x);

Interpol = zeros(nPlot,N);

for iPlot = 1:nPlot
    
    [isElement,pos] = ismember(x,xPlot(iPlot));
    
    if(any(isElement))
        
        pos = find(pos,1,'first');
        Interpol(iPlot,pos) = 1;
        
    else
        
        xPrevPos = find(x-xPlot(iPlot)<0,1,'last');
        xNextPos = find(x-xPlot(iPlot)>0,1,'first');
        
        xPrev = x(xPrevPos);
        xNext = x(xNextPos);
              
        alpha = (xPlot(iPlot) - xPrev)/(xNext - xPrev);
        
        Interpol(iPlot,xPrevPos) = 1 - alpha;
        Interpol(iPlot,xNextPos) = alpha;
        
    end
    
end


end