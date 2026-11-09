classdef WedgeCutSide < Polar_M1SpectralSpectral
    % A wedge of radius R_out cut off by a horizontal wall a distance h
    % below the Origin.  The angular range is such that the bottom corners
    % of the WedgeCut are where the circle of radii R_in and R_out centred
    % at Origin meet the wall.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R_out (or R) - outer radius of the wedge
    %  (R_in = 0) - inner radius of the wedge
    %  h - perpendicular distance from the Origin to the wall
    %  (Origin = [0,0]) - origin of the wedge


    properties
        R_in = 0;  % inner radius
        R_out      % outer radius
        h          % distance below origin
        leftRight = 'right'  % which side
        th1, th2   % construction angles
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        
        function this = WedgeCutSide(Geometry)
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            if(isfield(Geometry,'R'))
                this.R_out    = Geometry.R;
            else
                this.R_out    = Geometry.R_out;
            end
            if(isfield(Geometry,'R_in'))
                this.R_in     = Geometry.R_in;
            end

            this.h        = Geometry.h;
            this.leftRight = Geometry.leftRight;

            % invalid case
            if(this.h > this.R_in)
                exc = MException('Segment:Constructor','not yet implemented');
                throw(exc);
            end

            InitializationPts(this);

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            % map is not invertible at the corner
            this.invertibleMap = false; 

        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)
            int = ComputeIntegrationVector@Polar_M1SpectralSpectral(this);
            %Check Accuracy
            h     = this.h; %#ok<*PROP>
            Rin   = this.R_in;
            Rout  = this.R_out;
            thO   = 2*acos(h/Rout);
            thI   = 2*acos(h/Rin);
            area  = 0.5*(h*sqrt(Rout^2-h^2) - (Rin^2*thO/2 - Rin^2/2*(thI-sin(thI))));

            % display the error in the area if requested
            if(nargin == 2 && displayArea)
                if(area == 0)
                    disp(['WedgeCut: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['WedgeCut: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end
            end
        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)

            % intersection angles
            thO     = acos(this.h/this.R_out);
            thI     = acos(this.h/this.R_in);

            % adjust depending on which side we're on
            if(strcmp(this.leftRight,'right'))
                this.th1 =  3/2*pi+thI;
                this.th2 =  3/2*pi+thO;
            else
                this.th1 =  3/2*pi-thO;
                this.th2 =  3/2*pi-thI;
            end

            % construct points
            y2_kv   = LinearMap(x2,this.th1,this.th2);
            r_out   = (this.h)./abs(sin(y2_kv));
            y1_kv   = LinearMap(x1,this.R_in*ones(size(x1)),r_out);

            % Jacobian and Hessians
            if(nargout >= 3)
                J        = zeros(this.M,2,2);
                J(:,1,1) = (r_out-this.R_in)/2;
                J(:,2,2) = (this.th2 - this.th1)/2;
            end

            if(nargout >= 4)
                dH1        = zeros(this.M,2,2);
            end

            if(nargout >= 4)
                dH2        = zeros(this.M,2,2);
            end

        end

        function [x1,x2] = CompSpace(this,y1,y2) %#ok
            exc = MException('WedgeCutSide:CompSpace','not yet implemented');
            throw(exc);
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
            R_out = this.R_out;
            h = this.h;
            th1 = this.th1;
            th2 = this.th2;
            O1 = this.Origin(1);
            O2 = this.Origin(2);
            right = strcmp(this.leftRight,'right');

            % set sensible axis limits
            xLim = [-1.1 * R_out + O1, 1.1 * R_out + O1];
            yLim = [-1.1 * R_out + O2, 1.1 * R_out + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',15,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicators for the two angles
            tth1 = linspace(0, mod(th1,2*pi), 100);
            tth2 = linspace(0, mod(th2,2*pi), 100);
            R1 = R_in + (R_out-R_in)/4;
            R2 = R_in + (R_out-R_in)/3;
            plot(R1 * cos(tth1) + O1, R1 * sin(tth1) + O2, ...
                '-m', 'DisplayName', 'th1, th2');
            plot(R2 * cos(tth2) + O1, R2 * sin(tth2) + O2, ...
                '-m', 'HandleVisibility','off');
            % line at theta = 0
            plot([O1,xLim(2)],[O2,O2],':m','HandleVisibility','off');

            % plot inner and outer circles
            % R_in
            if(R_in > 0)
                tCircle = linspace(0, 2*pi, 100);
                plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                    '--k', 'HandleVisibility','off');
                p1 = [O1, O2];
                theta = 3*pi/4;
                p2 = [O1 + R_in * cos(theta), O2 + R_in * sin(theta)];
                dp = p2 - p1;
                plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
                text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2 + 0.02 * R_out, 'R in', ...
                    'FontSize',20,'Color','b')
            end

            % R_out
            tCircle = linspace(0, 2*pi, 100);
            plot(R_out*cos(tCircle) + O1, R_out*sin(tCircle) + O2, ...
                    '--k', 'DisplayName','Circles and wall');
            p1 = [O1, O2];
            theta = pi/4;
            p2 = [O1 + R_out * cos(theta), O2 + R_out * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2 - 0.02 * R_out, 'R out', ...
                'FontSize',20,'Color','b')
            
            % wall
            plot(xLim,[O2 - h, O2 - h],'--k','HandleVisibility','off');

            % plot indicator for h
            p1 = [O1, O2];
            p2 = [O1, O2 - h];
            dp = p2 - p1;
            plot([p1(1), p2(1)], [p1(2), p2(2)], '-b', 'DisplayName','h');
            text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2, ['h = ' num2str(h)], ...
                'FontSize',20,'Color','b')


            % plot the inner circle
            if(R_in > 0)
                tCircle = linspace(th1, th2, 100);
                plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                    '-r', 'HandleVisibility','off');
            end
            % plot the lines
            if(right)
                plot([R_in*cos(th2) + O1, R_out*cos(th2) + O1],[R_in*sin(th2) + O2, R_out*sin(th2) + O2], ...
                        '-r', 'DisplayName','Boundary');
                plot([R_in*cos(th1) + O1, R_out*cos(th2) + O1],[R_out*sin(th2) + O2, R_out*sin(th2) + O2], ...
                        '-r', 'HandleVisibility','off');
                plot([R_in*cos(th1) + O1, R_out*cos(th1) + O1],[R_in*sin(th1) + O2, R_out*sin(th1) + O2], ...
                        ':k', 'HandleVisibility','off');
            else
                plot([R_in*cos(th1) + O1, R_out*cos(th1) + O1],[R_in*sin(th1) + O2, R_out*sin(th1) + O2], ...
                        '-r', 'DisplayName','Boundary');
                plot([R_out*cos(th1) + O1, R_in*cos(th2) + O1],[R_in*sin(th2) + O2, R_in*sin(th2) + O2], ...
                        '-r', 'HandleVisibility','off');                
                plot([R_in*cos(th2) + O1, R_out*cos(th2) + O1],[R_in*sin(th2) + O2, R_out*sin(th2) + O2], ...
                        ':k', 'HandleVisibility','off');
            end

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
            geom.R_in = 2;
            geom.R_out = 3;
            geom.h = 1.9;
            geom.Origin = [-0.2,0.5];
            geom.leftRight = 'right';

        end

    end


end