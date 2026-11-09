classdef InfWedge < Polar_SpectralSpectral
    % An infinite wedge of a disc, between th1 and th2, centred at Origin.
    % Can also be used to create a wedge of an annulus by specifying R_in,
    % the inner radius.
    %
    % Inputs:
    %   N - [N1,N2] number of points in the two coordinates
    %   LR - mapping parameter for radial coordinate; half the points are
    %        in [R_in,R_in+LR]
    %   th1, th2 - the start and end angles
    %   (Origin = [0,0]) - origin of the circle

    properties
        R_in = 0 % inner radius
        LR       % mapping parameter
        th1,th2  % angle range
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = InfWedge(Geometry)
            this@Polar_SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            if(isfield(Geometry,'R_in'))
                this.R_in     = Geometry.R_in;
            end

            this.LR       = Geometry.LR;
            
            this.th1      = Geometry.th1;
            this.th2      = Geometry.th2;

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % set the origin
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)

            % use the standard integration vector
            int  = ComputeIntegrationVector@Polar_SpectralSpectral(this);

            area = Inf;

            % display the error in the area if requested
            if(nargin == 2 && displayArea)
                disp('InfWedge: has infinite area');
            end
        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % quotient map from [-1,1] to [R_in,inf]
            [r,dr,dx,ddx,dddx,ddddx] = QuotientMap(x1,this.LR,this.R_in);
        end
        function x1 = CompSpace1(this,r)
            % inverse quotient map
            x1  = InvQuotientMap(r,this.LR,this.R_in);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % linear map from [-1,1] to [th1,th2]
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.th1,this.th2);
        end
        function x2 = CompSpace2(this,y2)
            % inverse linear map
            x2 = InvLinearMap(y2,this.th1,this.th2);
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
            R_in = this.R_in; %#ok<*PROPLC>
            LR = this.LR;
            th1 = this.th1;
            th2 = this.th2;
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-3.1 * LR + O1, 3.1 * LR + O1];
            yLim = [-3.1 * LR + O2, 3.1 * LR + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circles and edges
            tCircle = linspace(th1, th2, 100);
            hold on
            if(R_in > 0)
                plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                    '-r', 'DisplayName', 'Boundary');
            end
            plot([R_in*cos(th1) + O1, 4*LR*cos(th1) + O1],[R_in*sin(th1) + O2, 4*LR*sin(th1) + O2], ...
                    '-r', 'HandleVisibility','off');
            plot([R_in*cos(th2) + O1, 4*LR*cos(th2) + O1],[R_in*sin(th2) + O2, 4*LR*sin(th2) + O2], ...
                    '-r', 'HandleVisibility','off');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',15,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicator for LR
            plot((LR+R_in)*cos(tCircle) + O1, (LR+R_in)*sin(tCircle) + O2, ...
                '--m', 'DisplayName', 'LR + R in');

            % arrow for R_in
            if(R_in > 0)
                theta = (th1+th2)/2;
                p1 = [O1, O2];
                p2 = [O1 + R_in * cos(theta), O2 + R_in * sin(theta)];
                dp = p2 - p1;
                plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
                text((p1(1) + dp(1)/2) + 0.02 * 3*LR, p1(2) + dp(2)/2 + 0.02 * 3*LR, 'R', ...
                    'FontSize',20,'Color','b')
            end

            % plot indicators for the two angles
            tth1 = linspace(0, mod(th1,2*pi), 100);
            tth2 = linspace(0, mod(th2,2*pi), 100);
            R1 = R_in + 0.2*LR;
            R2 = R_in + 0.3*LR;
            plot(R1 * cos(tth1) + O1, R1 * sin(tth1) + O2, ...
                '-m', 'DisplayName', 'th1, th2');
            plot(R2 * cos(tth2) + O1, R2 * sin(tth2) + O2, ...
                '-m', 'HandleVisibility','off');
            % line at theta = 0
            plot([O1,xLim(2)],[O2,O2],':m','HandleVisibility','off')

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

            geom.N = [10;12];
            geom.R_in = 4;
            geom.LR = 3;
            geom.th1 = pi/6;
            geom.th2 = 1.5*pi;
            geom.Origin = [-0.2,0.5];

        end

    end


end