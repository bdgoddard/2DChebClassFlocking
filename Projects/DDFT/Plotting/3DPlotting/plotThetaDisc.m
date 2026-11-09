function hFull = plotThetaDisc(theta,f,maxVal,minVal,x0,y0,r,cMap)
    N = length(theta);
    dTheta = pi/(N-1); % half width
    
    fRange = maxVal - minVal;
   
    %fVal = round( (f - minVal)/fRange * (N-1) ) + 1;
    
    m = 10;
    
    fVal = min(m,round((m-1)*(f-minVal)/fRange)+1);
    
    %cMap = colormap(gray(N));
    
    hFull = zeros(N,1);
    
    for iN = 1:N
        theta1 = theta(iN) - dTheta;
        theta2 = theta(iN) + dTheta;
        h = plotArc(theta1,theta2,x0,y0,r);        
        colour = cMap(fVal(iN),:);
        set(h,'facecolor',colour,'edgecolor','none');
        hold on
        hFull(iN) = h;
    end
   
    
        
        

end