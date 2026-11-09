classdef PeriodicHalfSpace < SpectralFourier
    % A half space which is periodic in the second cooodinate, with
    % length determined by y2Min and y2Max.
    % Alternatively, specify the length of the second side, L2.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  (y1Min = 0)   - left edge
    %  (y2Min = 0)   - bottom
    %  y2Max (or L2) - top
    %  L1            - mappring parameter for the first coordinate


    properties
        y1Min = 0  % left
        y2Min      % bottom (first periodic boundary)
        y2Max      % top (second periodic boundary_)
        L1 % mapping parameter for half infinite domain
        L2 % length of periodic domain
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = PeriodicHalfSpace(Geometry)
            this@SpectralFourier(Geometry.N(1),Geometry.N(2));

            % left
            this.y1Min = Geometry.y1Min;
            % bottom and top
            this.y2Min = Geometry.y2Min;
            this.y2Max = Geometry.y2Max;

            % set the mapping parameter for the infinite dimension
            this.L1 = Geometry.L1;

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

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % quotient map from [-1,1] to [y1Min,inf]
            [y1,dy1,dx,ddx,dddx,ddddx] = QuotientMap(x1,this.L1,this.y1Min,inf);
        end

        function x1 = CompSpace1(this,y1)
            % inverse quotient map
            x1  = InvQuotientMap(y1,this.L1,this.y1Min,inf);
        end

        function [y2,dy2,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % linear map from [0,1] to [y2Min,y2Max]
            % note that this is Fourier so defined on [0,1]
            [y2,dy2,dx,ddx,dddx,ddddx] = LinearMap01(x2,this.y2Min,this.y2Max);
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
            y1Min = this.y1Min; %#ok<*PROPLC>
            y2Min = this.y2Min;
            y2Max = this.y2Max;
            L1 = this.L1;
            L2 = y2Max - y2Min;

            % set sensible axis limits
            xLim = [y1Min - 0.1*L1, 3*L1];
            yLim = [y2Min - 0.1*L2, y2Max + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;

            hold on
            % plot the half space line
            plot([y1Min,y1Min],[y2Min,y2Max],'-r','DisplayName', 'Half space boundary');

            % plot the periodic lines
            plot([y1Min,xLim(2)],[y2Min,y2Min],'--r','DisplayName', 'Periodic boundary');
            plot([y1Min,xLim(2)],[y2Max,y2Max],'--r','HandleVisibility','off');

            % plot indicatorsfor L1
            plot([y1Min + L1, y1Min + L1],[y2Min,y2Max],'--m','DisplayName', 'Origin(1) + L1');

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
            geom.y1Min = -3;
            geom.L1 = 4;
            geom.y2Min = 0;
            geom.y2Max = 10;

        end

    end


end