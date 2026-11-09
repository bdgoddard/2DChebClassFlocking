classdef InfAnnulus < Polar_SpectralFourierNoOrigin
    % An infinite annulus with inner radius RMin and outer radius inf, centred
    % at Origin.
    % Note that this is defined in polar coordinates with y1 the radius (on
    % [RMin,inf]) and y2 the angle (on [0,2pi)).
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  RMin - the inner radius of the annulus
    %  L - a mapping parameter such that half of the radial points lie in
    %  [0,L]
    %  (Origin = [0,0]) - [O1,O2] the origin of the annulus

    properties
        RMin  % inner radius
        L     % mapping parameter
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        
        function this = InfAnnulus(Geometry)

            this@Polar_SpectralFourierNoOrigin(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.RMin = Geometry.RMin;

            if(isfield(Geometry,'L'))
                this.L    = Geometry.L;
            elseif(isfield(Geometry,'f'))
                this.L = FindOptimalL(@RadMap,this.Pts.x1,@Geometry.f,10);
            else
                error('InfAnnulus: No L or f given in input.');
            end

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % function if we want to ty to find the optimal L
            function r = RadMap(x,L)
                r = QuotientMap(x,L,this.RMin,inf);
            end
        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this)

            % standard integration vector
            int  = ComputeIntegrationVector@Polar_SpectralFourierNoOrigin(this);
            
            % provide the exact area
            area = Inf;
        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            % quotient map from [-1,1] to [RMin,inf]
            % such that half the points lie in [RMin,RMin + L]
            [r,dr,dx,ddx,dddx,ddddx] = QuotientMap(x,this.L,this.RMin,inf);
        end

        function xf = CompSpace1(this,r)
            % inverse quotient map
            xf = InvQuotientMap(r,this.L,this.RMin,inf);
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
            RMin = this.RMin; %#ok<*PROPLC>
            L = this.L; 
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * 3*(L+RMin) + O1, 1.1 * 3*(L+RMin) + O1];
            yLim = [-1.1 * 3*(L+RMin) + O2, 1.1 * 3*(L+RMin) + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot the circle at L+RMin
            tCircle = linspace(0, 2*pi, 100);
            plot((L+RMin)*cos(tCircle) + O1, (L+RMin)*sin(tCircle) + O2, ...
                '--m', 'DisplayName', 'Circle of radius L+R');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','Origin')

            % plot the boundary
            plot(RMin*cos(tCircle) + O1, RMin*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');

            % fix the axes
            axis equal
            box on
            set(gca,'XTick',[])
            set(gca,'YTick',[])            
            xlim(xLim);
            ylim(yLim);
            legend('Location','northwest');
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
            geom.RMin = 1;
            geom.L = 3;
            geom.Origin = [0.2;0.7];

        end

    end
    

end