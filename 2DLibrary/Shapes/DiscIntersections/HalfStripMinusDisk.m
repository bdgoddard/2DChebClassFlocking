classdef HalfStripMinusDisk < Polar_M1SpectralSpectral
    % Part of a half-infinite strip (either 'Left' or 'Right') cut off by a
    % disc of radius R centred at Origin on either the 'Top' or 'Bottom' edge.
    % The other edge of the strip is at y2Wall.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radius of the disc to be removed
    %  L1 - mapping parameter for radial coordinate; half the points lie in
    %       [R,R+L1]
    %  y2Wall - position of lower wall
    %  LeftRight - 'Left' or 'Right' for direction of strip
    %  TopBottom - 'Top' or 'Bottom' for whether strip is at
    %        the Top or Bottom of the disc.
    %  (Origin = [0,0]) - origin of the disc

    properties
        R  % radius of disc
        y2Wall  % location of slit wall
        L1  % mapping parameter
        LeftRight  % direction of slit - 'Left'/'Right'
        TopBottom  % location of Origin wall - 'Top'/'Bottom'
        th1, th2 % theta for DocumentationPlot
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = HalfStripMinusDisk(Geometry)

            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.y2Wall     = Geometry.y2Wall;
            this.LeftRight  = Geometry.LeftRight;

            if(isfield(Geometry,'Origin'))
                this.Origin     = Geometry.Origin;
            end
            if(isfield(Geometry,'TopBottom'))
                this.TopBottom = Geometry.TopBottom;
            end

            % ensure we're in polar coordinates
            this.polar = 'polar';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)

            % determine the angles of intersection of the disc and the
            % strip walls
            O    = ones(size(x1));
            d_th = asin(abs(this.Origin(2)-this.y2Wall)/this.R);

            if(strcmp(this.LeftRight,'Left'))
                if(strcmp(this.TopBottom,'Bottom'))
                    thMin = pi;
                    thMax = pi + d_th;
                else
                    thMin = pi - d_th;
                    thMax = pi;
                end
            elseif(strcmp(this.LeftRight,'Right'))
                if(strcmp(this.TopBottom,'Bottom'))
                    thMin = 2*pi - d_th;
                    thMax = 2*pi;
                else
                    thMin = 0;
                    thMax = d_th;
                end

            else
                error('HalfStripMinusDisk:Choose Left or Right.');
            end

            % store theta range in the shape
            this.th1 = thMin;
            this.th2 = thMax;

            % linear map in the angles
            [y2_kv,dy2] = LinearMap(x2,thMin,thMax);

            % rMax depends on theta
            rMax        = max(abs((this.Origin(2)-this.y2Wall)./sin(y2_kv)),this.R);

            % deal with infinity
            if(strcmp(this.LeftRight,'Left'))
                rMax(y2_kv == pi) = inf;
            else
                rMax((y2_kv == 2*pi) | (y2_kv == 0)) = inf;
            end
            % apply a quotient map
            rd                = rMax-this.R;
            L1_r              = this.L1*rd./(this.L1*3+rd);
            L1_r(rd==inf)     = this.L1;
            [y1_kv,dy1]       = QuotientMap(x1,L1_r,O*this.R,max(rMax,this.R));

            % Jacobian and Hessian
            n = this.N1*this.N2;
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = dy1;
                J(:,2,2) = dy2;
            end

            if(nargout >= 4)
                dH1        = zeros(n,2,2);
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);
            end
        end

        function CompSpace(~,~,~)
            % computational space not yet implemented
            exc = MException('Segment:CompSpace','not yet implemented');
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
            L1 = this.L1; %#ok<*PROPLC>
            R  = this.R;
            y2Wall = this.y2Wall;
            th1 = this.th1;
            th2 = this.th2;
            O1 = this.Origin(1);
            O2 = this.Origin(2);
            L2 = abs(O2 - y2Wall);
            right = strcmp(this.LeftRight,'Right');
            top = strcmp(this.TopBottom,'Top');

            % set sensible axis limits
            xLim = [-1.1 * 3*L1 + O1, 1.1 * 3*L1 + O1];
            if(top)
                yLim = [O2 - 0.5*L2, y2Wall + 0.5*L2];
            else
                yLim = [y2Wall - 0.5*L2, O2 + 0.5*L2];
            end

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator line for L1+R
            tCircle = linspace(0, 2*pi, 100);
            plot((L1+R)*cos(tCircle) + O1, (L1+R)*sin(tCircle) + O2, ...
                '--m', 'DisplayName','L1 + R');

            % plot circle
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName','Walls and Disc');
            p1 = [O1, O2];
            if(top)
                theta = 3*pi/4;
            else
                theta = -3*pi/4;
            end
            p2 = [O1 + R * cos(theta), O2 + R * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * 3*L1, p1(2) + dp(2)/2 + 0.02 * 3*L1, 'R', ...
                'FontSize',20,'Color','b')

            % plot walls
            plot(xLim,[O2,O2],'--k', 'HandleVisibility','off');
            plot(xLim,[y2Wall,y2Wall],'--k', 'HandleVisibility','off');

            % plot boundary

            if(right)
                plot([R + O1,xLim(2)],[O2,O2],'-r', 'DisplayName', 'Boundary');
                plot([O1 + R*cos(th1),xLim(2)],[y2Wall,y2Wall],'-r', 'HandleVisibility','off');
            else
                plot([xLim(1),O1 - R],[O2,O2],'-r', 'DisplayName','Boundary');
                plot([xLim(1),O1 + R*cos(th2)],[y2Wall,y2Wall],'-r', 'HandleVisibility','off');
            end

            % plot curved edge
            tCircle = linspace(th1, th2, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'HandleVisibility','off');

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

            geom.N = [20;12];
            geom.R = 1;
            % geom.y2Wall = -0.6;
            % geom.TopBottom = 'Bottom';
            geom.L1 = 1;
            geom.y2Wall = 0.6;
            geom.TopBottom = 'Top';
            geom.LeftRight = 'Left';
            geom.Origin = [0,0];

        end

    end


end