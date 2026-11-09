classdef PeriodicInfBox < SpectralFourier
    % A (half) periodic box which is infinite in the first coordinate.
    % The box is periodic in the second coordinate, with length determined by
    % y2Min and y2Max.
    % Alternatively, specify the length of the second side, L2.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  (y2Min = 0)   - bottom
    %  y2Max (or L2) - top
    %  L1            - mapping parameter for the first coordinate

    properties
        y2Min = 0 % bottom (first periodic boundary)
        y2Max     % top (second periodic boundary_)
        L1 % note that this isn't the size of the y1 domain; that's infinite,
           % rather, we have that half the points lie in [-L,L]
        L2 % length of periodic dimension
    end
    
    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = PeriodicInfBox(Geometry)
            this@SpectralFourier(Geometry.N(1),Geometry.N(2));

            % set the mapping parameter for the infinite dimension
            this.L1 = Geometry.L1;

            % set top, and bottom
            this.y2Min = Geometry.y2Min;
            this.y2Max = Geometry.y2Max;

            % determine the periodic length
            this.L2 = this.y2Max - this.y2Min;

            % set that we're in Cartesian coordinates
            this.polar = 'cart';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------
        
        function int = ComputeIntegrationVector(this)
            % overridden from SpectralFourier to deal with infinite entries in
            % the integration vector

            int = ComputeIntegrationVector@SpectralFourier(this);
            int(~isfinite(int)) = 0; % deal with points at \pm \infty
            this.Int = int;
        end


        function M_IntConv = ComputeIntegrationConvolutionMatrix(this,fInt,fConv,shapeParams)
            % Computes the matrix that does a weighted integral in the
            % first coordinate (with function fInt) and a convolution in
            % the second coordinate (with fConv).  The convolution can be
            % performed on a subdomain specified by shapeParams.

            if(nargin<4)
                shapeParams = struct;
            end

            % Weighted integration in y1 (spectral)

            % set up spectral line
            geom1.N = this.N1;
            geom1.L = this.L1;
            sLine = InfSpectralLine(geom1);

            % compute integration weights
            Int1 = sLine.ComputeIntegrationVector;
            ws = sLine.Pts.y;

            gLine = fInt(ws);
            %TOM: removed non-finite points (11/08/20)
            gLine(~isfinite(gLine)) = 0;
            I1Line = Int1.* gLine';
            I1Rep = repmat(I1Line,[this.N1,1]);
            I1 = sparse(kron(I1Rep,eye(this.N2)));

            % Convolution in y2 (Fourier)

            % set up Fourier line
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);

            % compute convolution (which requires the integration vector)
            fLine.ComputeIntegrationVector;
            % format: C2 = fLine.ComputeConvolutionMatrix(this,f,shapeParams,parent)
            C2Line = fLine.ComputeConvolutionMatrix(fConv,shapeParams);

            C2 = sparse(kron(eye(this.N1),C2Line));

            % Full integral
            M_IntConv = C2*I1;
        end

        function [Int1,Int2] = Compute1DIntegrationMatrices(this)
            % Partial integration in each of the two coordinates

            % Spectral integration
            geom1.N = this.N1;
            geom1.L = this.L1;
            sLine = InfSpectralLine(geom1);
            I1 = sLine.ComputeIntegrationVector;
            Int1 = sparse(kron(I1,eye(this.N2)));

            % Fourier integration
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);
            I2 = fLine.ComputeIntegrationVector;
            Int2 = sparse(kron(eye(this.N1),I2));
        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % square root map from [-1,1] to [-inf,inf]
            [y,dy,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);
        end

        function xf = CompSpace1(this,y1)
            % inverse square root map
            xf  = InvSqrtMap(y1,this.L1,inf);
        end

        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % linear map from [0,1] to [y2Min,y2Max]
            % note that this is Fourier so defined on [0,1]
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(x2,this.y2Min,this.y2Max);
        end

        function x2 = CompSpace2(this,y2)
            % inverse linear map
            x2 = InvLinearMap01(y2,this.y2Min,this.y2Max);
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
            y2Min = this.y2Min; %#ok<*PROPLC>
            y2Max = this.y2Max;
            L1 = this.L1;
            L2 = y2Max - y2Min;

            % set sensible axis limits
            xLim = [-3*L1, 3*L1];
            yLim = [y2Min - 0.1*L2, y2Max + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;

            % plot the box
            hold on
            plot(xLim,[y2Min,y2Min],'--r','DisplayName', 'Periodic boundary');
            plot(xLim,[y2Max,y2Max],'--r','HandleVisibility','off');

            % plot indicators for L1
            plot([-L1,-L1],[y2Min,y2Max],'--m','DisplayName', 'L1');
            plot([L1,L1],[y2Min,y2Max],'--m','HandleVisibility','off');

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

            geom.N = [20;12];
            geom.L1 = 4;
            geom.y2Min = 0;
            geom.y2Max = 10;

        end

    end

end