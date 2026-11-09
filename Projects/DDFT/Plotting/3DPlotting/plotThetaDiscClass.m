function plotThetaDiscClass(f,x0,y0,r,m)
    
    % number of spatial points to plot at
    NPts = size(f,2);
    
    % set up annulus with origin [0;0]
    Nt = size(f,1);
    Nr = 10;
    geom.N = [Nr;Nt];
    geom.R_in = r/4;
    geom.R_out = r;

    interp1 = (-1:0.1:1)';
    interp2 = (0:0.01:1)';
    
    annulus = Annulus(geom);
    
    % compute interpolation matrix
    Interp = ComputeInterpolationMatrix(annulus,interp1,interp2,true,true);
    InterPol = Interp.InterPol;
    

    % determine bounds on f - used to set colour range
    cMin = 0;
    cMax = 0;
    
    for iPt = 1:NPts
        fPt = f(:,iPt);
        fPt = kron(ones(Nr,1),fPt);
        fPtInterp = InterPol*fPt;
        ptMin = min(min(fPtInterp));
        ptMax = max(max(fPtInterp));
        cMin = min(cMin,ptMin);
        cMax = max(cMax,ptMax);
    end

    options = 'color';
    optDetails.labels = false;

    for iPt = 1:NPts
        % map f onto [0,m]
        fPt = f(:,iPt);
        fPt = min(m,(m-1)*(fPt-cMin)/(cMax-cMin) +1); 
        
        % replicate over radial coordinate
        fPt = kron(ones(Nr,1),fPt);
        x0Pt = x0(iPt);
        y0Pt = y0(iPt);
        
        % shift origin to plotting point
        annulus.Origin = [x0Pt;y0Pt];
        
        % plot data
        annulus.plot(fPt,options,optDetails);
        hold on
    end
        

end