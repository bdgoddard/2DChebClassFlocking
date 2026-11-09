classdef Annulus < Polar_SpectralFourierNoOrigin
    % An annulus with inner radius RMin and outer radius RMax, centred
    % at Origin.
    % Note that this is defined in polar coordinates with y1 the radius (on
    % [RMin,RMax]) and y2 the angle (on [0,2pi)).
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  RMin (or R_in)   - the inner radius of the annulus
    %  RMax (or R_out)  - the inner radius of the annulus
    %  (Origin = [0,0]) - [O1,O2] the origin of the annulus

    properties
        RMin  % inner radius
        RMax  % outer radius
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        
        function this = Annulus(Geometry)

            this@Polar_SpectralFourierNoOrigin(Geometry.N(1),Geometry.N(2));

            % set parameters
            if(isfield(Geometry,'R_in'))
                this.RMin   = Geometry.R_in;
                this.RMax   = Geometry.R_out;
            else
                this.RMin   = Geometry.RMin;
                this.RMax   = Geometry.RMax;
            end

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

        function [int,area] = ComputeIntegrationVector(this)

            % standard integration vector
            int = ComputeIntegrationVector@Polar_SpectralFourierNoOrigin(this);

            % provide exact area
            area = pi*(this.RMax^2-this.RMin^2);
        end


        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            % standard linear map from [-1,1] to [-R,R]
            [r,dr,dx,ddx,dddx,ddddx] = LinearMap(x,this.RMin,this.RMax);
        end

        function xf = CompSpace1(this,r)
            % inverse linear map
            xf = InvLinearMap(r,this.RMin,this.RMax);
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
        % Plotting
        %------------------------------------------------------------------

        function PlotBorders(this)
            % plot the borders of the annulus
            % complements PlotBound of Shape
            
            % get the indices
            ComputeIndices(this);

            % find the inner circle and plot it
            rIn = this.Pts.y1_kv(this.Ind.inR);
            geomIn.R = rIn(1);
            geomIn.N = this.Pts.N2;
            cIn = Circle(geomIn);

            % find the outer circle and plot it
            rOut = this.Pts.y1_kv(this.Ind.outR);
            geomOut.R = rOut(1);
            geomOut.N = this.Pts.N2;
            cOut = Circle(geomOut);

            cIn.PlotFunction;
            cOut.PlotFunction;

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
            RMax = this.RMax;
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * RMax + O1, 1.1 * RMax + O1];
            yLim = [-1.1 * RMax + O2, 1.1 * RMax + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circles
            tCircle = linspace(0, 2*pi, 100);
            hold on
            plot(RMin*cos(tCircle) + O1, RMin*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            plot(RMax*cos(tCircle) + O1, RMax*sin(tCircle) + O2, ...
                '-r', 'HandleVisibility','off');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicator arrows
            % arrow for RMin
            p1 = [O1, O2];
            theta = pi/6;
            p2 = [O1 + RMin * cos(theta), O2 + RMin * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * RMax, p1(2) + dp(2)/2 - 0.02 * RMax, 'RMin', ...
                'FontSize',20,'Color','b')
            % arrow for RMax
            theta = 2*pi/6;
            p2 = [O1 + RMax * cos(theta), O2 + RMax * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * RMax, p1(2) + dp(2)/2 - + 0.02 * RMax, 'Rmax', ...
                'FontSize',20,'Color','b')

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
            geom.RMin = 1;
            geom.RMax = 3;
            geom.Origin = [0.2;0.7];

        end

    end

end