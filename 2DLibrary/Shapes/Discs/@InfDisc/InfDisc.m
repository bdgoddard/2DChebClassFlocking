classdef InfDisc < Polar_SpectralFourier
    % A disc of infinite radius, centred at Origin.  I.e., the plane in
    % polar coordinates.
    % Note that this is defined in polar coordinates with y1 the radius (on
    % (0,inf]) and y2 the angle (on [0,2pi)).
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L - a mapping parameter such that half of the radial points lie in
    %  [0,L]
    %  (Origin = [0,0]) - [O1,O2] the origin of the disc

    properties
        L  % mapping parameter
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = InfDisc(Geometry)

            this@Polar_SpectralFourier(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L = Geometry.L;

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
        end
    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams,parent)

            % force us to use the standard convolution
            if((nargin == 4) && parent)
                M_conv = ComputeConvolutionMatrix@Polar_SpectralFourier(this,f);
                return;
            end

            %{
            Strategy:
            We wish to compute \int_D(x) f(x-y) g(y) dy, for some domain D
            (e.g. disc, annulus) centred at x.
            Let z = y-x then dt = dx, y = x + z and 
            y \in D(x) => z + x \in D(x) or z \in D(0) and we get
            \int_D(0) f(-z) g(x+z) dz
            
            In matlab this is
            M_conv = Int_D(0) .* f(-z) .* IP_D(x)
            %}

            disp('Computing Convolution matrices...');

            % shorthand some parameters
            n1 = this.N1;
            n2 = this.N2;
            y1_kv = this.Pts.y1_kv;
            y2_kv = this.Pts.y2_kv;

            % determine which subshape we're using
            if(isfield(shapeParams,'R'))
                % disc
                subShape = Disc(shapeParams);
                rRange = (-1:0.05:1)';
            elseif(isfield(shapeParams,'RMin'))
                % annulus
                if(isfield(shapeParams,'RMax'))
                    subShape = Annulus(shapeParams);
                    rRange = (-1:0.05:1)';
                else
                    subShape = InfAnnulus(shapeParams);
                    rRange = (-1:0.05:0.9)';
                end
            else
                % default to infinite disc
                subShape = InfDisc(shapeParams);
                rRange = (-0.9:0.05:0.9)';
            end

            % Integration and interpolation on D(0)
            Int_i    = subShape.ComputeIntegrationVector();
            Interp   = subShape.ComputeInterpolationMatrix(rRange,(0:0.02:1)',true,false);

            % f(-z)
            invSubShapePts = invertPts(subShape.Pts,'polar');
            fP       = f(GetDistance(this,invSubShapePts.y1_kv,invSubShapePts.y2_kv));

            % deal with multiple kernels
            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end
            nf = size(fP,2);

            % replicate the integration vector
            IntRep = Int_i(ones(1,nf),:);

            % preallocate the convolution matrix
            M_conv = zeros(n1*n2,n1*n2,nf);

            % plot the kernels
            figure
            subShape.Interp=Interp;
            subShape.plot(fP,'SC');

            % compute the convolution matrices
            hw = waitbar(0,'Computing Convolution matrices...');

            % loop over all points in the InfDisc
            for i=1:(n1*n2)

                % progress update
                waitbar(i/(n1*n2),hw);

                % the current point
                y10 = y1_kv(i);
                y20 = y2_kv(i);

                % compute points in D(x) = (y10,y20))
                subConvShapePts = shiftPointsPolar(subShape.Pts,y10,y20,'polar');
                % and the interpolation onto D(x) from the full domain
                IP = SubShapePts(this,subConvShapePts);

                % compute the convolution 
                IntFTemp = (IntRep.*fP.')*IP;
                M_conv(i,:,:) = IntFTemp.';

            end

            % close the weightbar
            close(hw);

            % deal with any issues caused by nans
            M_conv(isnan(M_conv)) = 0;

            % reshape so the output is of the correct form
            if(nDim>2)
                nfTemp = size(fP_old);
                nfTemp = nfTemp(2:end);
                M_conv = reshape(M_conv,[ n1*n2, n1*n2, nfTemp ]);
            end


        end

        % test of convolution methodology against known solutions
        [M_conv,M_conv_New] = ComputeConvolutionMatrix_Test(this,shapeParams)

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            % square root map from [-1,1] to [-inf,inf]
            % such that half the points lie in [0,L]
            % note that Polar_SpectralFourier will cut this down to only
            % the positive radius values, after having doubled N1
            % so the points are singly defined, but there is more effort
            % needed to compute the differentiation and integration           
            [r,dr,dx,ddx,dddx,ddddx] = SqrtMap(x,this.L,inf);
        end

        function xf = CompSpace1(this,r)
            % inverse inverse square root map
            xf = InvSqrtMap(r,this.L,inf);
        end

        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(~,xT)
            % linear map from [0,1] to [0,2pi]
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end

        function xf = CompSpace2(~,th)
            % inverse linear map
            xf = InvLinearMap01(th,0,2*pi);
        end


        %------------------------------------------------------------------
        % Documentation plotting
        %------------------------------------------------------------------

        function [hf,ha,plotFile] = DocumentationPlot(this,opts)

            if(nargin<2)
                opts = struct;
            end

            if(isfield(opts,'savePlot'))
                savePlot = opts.savePlot;
            else
                savePlot = false;
            end

            % shorthand parameters
            L = this.L; %#ok<*PROPLC>
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * 3*L + O1, 1.1 * 3*L + O1];
            yLim = [-1.1 * 3*L + O2, 1.1 * 3*L + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circle at L
            tCircle = linspace(0, 2*pi, 100);
            hold on
            plot(L*cos(tCircle) + O1, L*sin(tCircle) + O2, ...
                '--m', 'DisplayName', 'Circle of radius L');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','Origin')

            % fix the axes
            axis equal
            box on
            set(gca,'XTick',[])
            set(gca,'YTick',[])            
            xlim(xLim);
            ylim(yLim);
            legend('Location','best');
            title(mfilename);

            % save image to file if required
            if(savePlot)
                if(isfield(opts,'plotFile'))
                    plotFile = opts.plotFile;
                else
                    plotFile = [mfilename '.png'];
                end

                saveDir = ['DocumentationImages' filesep];
                if(~exist(saveDir,'dir'))
                    mkdir(saveDir);
                end

                plotFile = [saveDir plotFile];
                save2png(plotFile,hf);
                close(hf);
            else
                plotFile = [];
            end
        end

    end

    methods (Static)

        function geom = DocumentationPlotDefaultGeometry()
            % set a sensible default geometry

            geom.N = [14;12];
            geom.L = 3;
            geom.Origin = [0.2,0.7];

        end

    end


end