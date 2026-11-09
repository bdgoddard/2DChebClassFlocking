function [fMin,fMax] = plotThetaAnnulus(f,x0,y0,r,m)
    
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
    fMin = 0;
    fMax = 0;

    for iPt = 1:NPts
        fPt = f(:,iPt);
        fPt = kron(ones(Nr,1),fPt);
        fPtInterp = InterPol*fPt;
        ptMin = min(min(fPtInterp));
        ptMax = max(max(fPtInterp));
        fMin = min(fMin,ptMin);
        fMax = max(fMax,ptMax);
    end

    % determine bounds on f - used to set colour range
    fMinc = 0;
    fMaxc = 0;

    geomc.R = 1;
    geomc.N = size(f,1);
    circle = Circle(geomc);
    
    Interpc = ComputeInterpolationMatrix(circle,interp2);
    
    InterPolc = Interpc.InterPol;
    
    for iPt = 1:NPts
        fPt = f(:,iPt);
        fPtInterp = InterPolc*fPt;
        ptMin = min(min(fPtInterp));
        ptMax = max(max(fPtInterp));
        fMinc = min(fMinc,ptMin);
        fMaxc = max(fMaxc,ptMax);
    end

    fMin
    fMinc
    
    fMax
    fMaxc
    
    
    options = 'color';
    optDetails.labels = false;

    for iPt = 1:NPts
        % map f onto [0,m]
        fPt = f(:,iPt);
        fPt = min(m,(m-1)*(fPt-fMin)/(fMax-fMin) +1); 
        
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