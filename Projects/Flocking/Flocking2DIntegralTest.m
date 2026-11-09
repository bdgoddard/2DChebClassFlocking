function data = Flocking2DIntegralTest(gamma)

%     y1 is Fourier - position
%     y2 is Fourier - angle
    
    disp('** Flocking2DIntegrationTest **');

    AddPaths;
   
    % set domain parameters

    % space
    xMin = 0;
    xMax = 1;
    % orientation
    Ltheta = 2*pi; % length of the domain	
    
    % set up the 2D domain and compute the required operators
    shape.y1Min = xMin;
    shape.y1Max = xMax;	
    N1          = 40;
    shape.y2Min = 0;
    shape.y2Max = Ltheta;
    N2      = 40;
    shape.N = [N1;N2];  
              
    aBox = BoxFourierFourier(shape);
    [Pts,Diff,Int,Ind] = aBox.ComputeAll();    

    % intepolation for plotting
    Interp             = aBox.ComputeInterpolationMatrix(...
                                     (0:0.025:1)',(0:0.025:1)',true,true);  

    % shorthands to make notation clearer
    x = Pts.y1_kv;
    theta = Pts.y2_kv;

    % parameters for the convolution
    R = 0.2;  % this is epsilon
    convParams1.N = 40; % number of points to use when computing the convolution integral on the interval
    convParams1.yMax = R;

    convParams2 = struct; % don't apply any cutoff in the orientation integral

    % compute the convolution integrals
    % M_ConvConv: the double convolution = C2*C1
    % C1: the spatial convolution
    % C2: the orientation convolution
    % C1 and C2 can be use to compute the 1D convolutions, e.g., for the
    % local / non-conservative case
    [M_ConvConv,C1,C2] = aBox.ComputeConvolutionConvolutionMatrix(@convFnSpace,@convFnOrientation,convParams1,convParams2);
    [M_ConvConvD,CD1,CD2] = aBox.ComputeConvolutionConvolutionMatrix(@convFnSpace,@convFnSpace,convParams1,convParams2);


    % % Test against an exact (double) convolution
    % alpha = 2*pi;
    % psi = sin(alpha * x + theta);
    % % angular integral with K = - cos(theta) gives 
    % % - pi sin(alpha*x + theta)
    % % spatial integral from a to b gives pi/alpha [ cos( alpha b + theta )
    % % - cos( alpha a + theta)
    % 
    % C = M_ConvConv*psi;
    % xm = x - R;
    % xp = x + R;
    % CExact = pi/alpha * ( cos(alpha * xp + theta) - cos(alpha * xm + theta) );
    % 
    % figure('Position',[100,100,1200,800])
    % subplot(1,3,1)
    % aBox.plot(C);
    % subplot(1,3,2)
    % aBox.plot(CExact);
    % subplot(1,3,3);
    % aBox.plot(C-CExact);
    % 
    % sqrt( Int*((C - CExact).^2) / (Int * (CExact).^2 ))
    % 
    % data = struct;
    % return


    function f = convFnSpace(x)
        % spatial integration weight (not taking into account the cutoff)
        f = ones(size(x));
    end

    function f = convFnOrientation(theta)
        % orientation integration weight (not taking into account the cutoff)
        f = sin(theta);
    end

    % differentiation matrices
    Dx = Diff.Dy1;
    Dtheta = Diff.Dy2;
    DDtheta = Diff.DDy2;
        
    % get initial condition and normalize
    psiInitial = initialCondition(x,theta);
    % Normalise
    psiInitial = psiInitial/(Int*psiInitial);
               
    function psi0 = initialCondition(x,theta)
        % this is arbtirary (note it needs to be periodic in both
        % dimensions
        psi0 = cos(theta) .* sin(2*pi*x) + 2;
        Z = Int*psi0;
        psi0= psi0/Z;
    end

   
    % Solve the PDE
    mM  = ones(size(x));
    opts = odeset('RelTol',10^-4,'AbsTol',10^-4,'Mass',diag(mM));

    % maximum time
    tMax = 250;

    % parameters from the model
    v0 = 0.1;
    if(nargin==0)
        gamma = 0;
    end       
    Gamma = 1;

    % do the PDE solve
    [outTimes,psi_t] = ode15s(@flockingPDE, [0:tMax/100:tMax],psiInitial,opts);

    % transpose so that the second dimension is time
    psi_t = psi_t';

    % make a movie
    figure
    

    title = ['Density $\rho(x, \theta)$. $R = ' num2str(R), '$, $\gamma = ' num2str(gamma) '$, $D(u) = -\cos(u)$, partial normalisation in $\theta$'];

    optDetails.xLabel = 'Space';
    optDetails.yLabel = 'Orientation';
    optDetails.fontsize = 26;  % set font size (20 is the standard)
    optDetails.title = title;  % set the title in the movie

    optDetails.pdfDir = 'TempMoviePDFs/'; % directory to store the pdfs (need the /)
    optDetails.doPdfs = false; % make pdfs of each frame

    optDetails.doMp4 = false; % make mp4 movie
    optDetails.mp4File = 'TestMovie.mp4'; % movie file name

    optDetails.zLim = [0,2];
    % optDetails.edgecolor = 'none';
    
    aBox.makeMovie(psi_t,{},optDetails);

    save2pdf('FinalFigure.pdf',gcf);
    
    % output the data
    data.psi = psi_t;
    data.aBox = aBox;


    % the PDE
    function dpsidt = flockingPDE(t,psi)
       

        %interaction = (M_ConvConv * psi); %./(M_ConvConvD * psi);

        numerator   = M_ConvConv * psi;    % spatial sine(θ-θ') convolution
        denominator = CD1       * psi;    % angular integral ρ(x) (no spatial averaging)

        %avoid division by zero
        denominator = denominator + 1e-12;

        interaction = numerator; %./denominator;

        dpsidt = - v0 * cos(theta) .* Dx * psi ...
                 + gamma * Dtheta * (psi .* (interaction)) ...
                 + Gamma * DDtheta * psi;
    end

end