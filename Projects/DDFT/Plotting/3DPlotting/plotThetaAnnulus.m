function plotThetaAnnulus(f,x0,y0,opts)
    
    r    = opts.r;
    fMin = opts.fMin;
    fMax = opts.fMax;
    m    = opts.m;

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

    options = 'color';
    optDetails.labels = false;

    for iPt = 1:NPts
        % map f onto [0,m]
        fPt = f(:,iPt);
        fPt = min(m,(m-1)*(fPt-fMin)/(fMax-fMin) +1); 
        
        % replicate over radial coordinate
        fPt = kron(ones(Nr,1),fPt);

        % shift origin to plotting point
        annulus.Origin = [x0(iPt);y0(iPt)];
        
        % plot data
        annulus.plot(fPt,options,optDetails);
        hold on
    end
        

end