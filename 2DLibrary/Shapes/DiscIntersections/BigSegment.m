classdef BigSegment < handle
    % A 'big' segment of a disc of radius R, created by a
    % horizontal/vertical chord at perpendicular position Wall_Y.
    % Note that this is the complement of Segement.
    %
    % This is formed from a triangle and wedge (the part of a disc between
    % two angles).
    %
    % Inputs:
    %   NT - [NT1,NT2] number of points for the triangle in the two coordinates
    %   NW - [NW1,NW2] number of points for the wedge in the two
    %   coordinates (radial and angular)
    %   R - radius of the disc
    %   Wall_VertHir - 'horizontal' or 'vertical' orientation of the wall
    %   Wall_Y - coordiate of the wall in the perpendicular direction
    %   Wall - one of 'N', 'S', 'E', 'W' denoting the compass direction of
    %   the wall from the origin.  Only used when we have exactly half a
    %   disc
    %   (Origin = [0,0]) - origin of the disc
    %   (sphere = false) - true implements 'spherical integration'


    properties
        Origin = [0;0]; % origin of the disc
        R % radius of disc
        Wall_Y % perpendicular position of wall wrt Origin
        W % a Wedge Shape
        T % a Triangle Shape
        h % distance of wall from Origin
        Int % integration vector on whole shape
        Pts % points
        sphere = false % whether to do spherically-weighted integration
        horizontal % horizontal/vertical flag
        invertibleMap = false; % map is not invertible due to Triangle
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = BigSegment(Geometry)

            % set parameters
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            if(isfield(Geometry,'sphere'))
                this.sphere    = Geometry.sphere;
                wedgeSh.sphere = Geometry.sphere;
            end

            this.R      = Geometry.R;
            this.Wall_Y = Geometry.Wall_Y;

            % determine th1 and th2 which describe the domain of the wedge
            if(isfield(Geometry,'Wall_Y'))
                % validate vertical/horizontal choice
                horizontal = strcmp(Geometry.Wall_VertHor,'horizontal');
                vertical   = strcmp(Geometry.Wall_VertHor,'vertical');
                if(~(horizontal || vertical))
                    exc = MException('SegmentWall:Constructor', 'Wall_VertHor must be horizontal or vertical');
                    throw(exc);
                end
                this.horizontal = horizontal;

                % determine the distance from the wall
                % and the two angles to the edges of the wedge
                if(horizontal)
                    h = Geometry.Wall_Y - Geometry.Origin(2);

                    if(h > this.R)
                        this.throwHRError();
                    end

                    Geometry.th1 = pi - asin(h/this.R);
                    Geometry.th2 = asin(h/this.R);
                else
                    h = Geometry.Wall_Y - Geometry.Origin(1);

                    if(h > this.R)
                        this.throwHRError();
                    end

                    Geometry.th1 = acos(h/this.R);
                    Geometry.th2 = -acos(h/this.R);
                end
            end

            % neaten up th1 and th2 to make the case analysis more
            % straightforward
            Geometry.th1 = mod(Geometry.th1,2*pi);
            Geometry.th2 = mod(Geometry.th2,2*pi);

            % ensure that the thetas reflect the geometry
            if(mod(Geometry.th1 - Geometry.th2,2*pi) < pi)
                % standard case
                th1 = Geometry.th1;
                th2 = Geometry.th2;

                % deal with the edge case where the wall coincides with the
                % origin
                % note: only need to specify Geometry.Wall in this case
            elseif(mod(Geometry.th1 - Geometry.th2,2*pi) == pi)

                if(strcmp(Geometry.Wall,'N') || strcmp(Geometry.Wall,'E'))
                    th1 = Geometry.th1;
                    th2 = Geometry.th2;
                end
                if(strcmp(Geometry.Wall,'S') || strcmp(Geometry.Wall,'W'))
                    th1 = Geometry.th2;
                    th2 = Geometry.th1;
                end

            else
                % need to redefine ordering in this case to give the
                % correctly oriented wedge

                th1 = Geometry.th2;
                th2 = Geometry.th1;

            end

            % set up the wedge
            wedgeSh.th1    = th1;
            wedgeSh.th2    = th2;
            wedgeSh.R_out  = this.R;
            wedgeSh.N      = Geometry.NW;
            this.W         = Wedge(wedgeSh);
            % include the Origin in the shape for plotting
            this.W.Origin  = this.Origin;
            % get the true points including the Origin (inherited from
            % Shape)
            wCartPts       = this.W.GetCartPts();

            % set up the triangle
            r   = Geometry.R;
            Y = [0,0;
                cos(th1)*r,sin(th1)*r;
                cos(th2)*r,sin(th2)*r];
            triangleSh.Y = Y;
            triangleSh.N = Geometry.NT;

            this.T = Triangle(triangleSh);
            % include the Origin in the shape for plotting
            this.T.Origin = this.Origin;
            % get the true points including the Origin (inherited from
            % Shape)
            tCartPts       = this.T.GetCartPts();

            % concatenate the points
            % this.Pts  = struct('y1_kv',[wCartPts.y1_kv;this.T.Pts.y1_kv],...
            %     'y2_kv',[wCartPts.y2_kv;this.T.Pts.y2_kv]);
            this.Pts  = struct('y1_kv',[wCartPts.y1_kv;tCartPts.y1_kv],...
                'y2_kv',[wCartPts.y2_kv;tCartPts.y2_kv]);

            % define h, the distance to the wall
            this.h = r*abs(cos((th1-th2)/2));
        end

        %------------------------------------------------------------------
        % Shape-like functions
        %------------------------------------------------------------------

        % Cartesian points
        function ptsCart = GetCartPts(this)
            ptsCart       = this.Pts;
        end

        function [int,area] = ComputeIntegrationVector(this,displayArea)

            % integration weights for the triangle
            [Tint,~] = this.T.ComputeIntegrationVector();

            % integration weights for the wedge
            if(this.sphere) % if we're spherically weighting

                % inputs are (t1Odd, t2Odd) - see Wedge
                [wint,~]  = this.W.ComputeIntegrationVector(true,false);

                % combine the integration weights
                int       = [wint,Tint];

                % do the spherical weighting
                y1s  = this.Pts.y1_kv;
                y2s  = this.Pts.y2_kv;
                int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';

                % update Int
                this.Int = int;

                % compute the exact area of the whole shape
                ht   = this.R - abs(this.h);
                th = acos(abs(this.h)/this.R);
                area1 = 4/3*(pi-th)*this.R^3;

                % and of the wedge and triangle
                intW  = int(1:length(this.W.Int));
                area = 4/3*pi*this.R^3 - pi*ht^2/3*(3*this.R - ht);
                area2 = area - area1;
                intT = int(length(this.W.Int)+1:end);

                % display the error in the individual areas if requested
                if(nargin == 2 && displayArea)
                    disp(['Error of integration of (spherical) Wedge (ratio):',num2str(1-abs(sum(intW)/area1))]);
                    disp(['Error of integration of (spherical) Wedge (abs):',num2str(abs(sum(intW)-area1))]);

                    disp(['Error of integration of (spherical) Triangle (ratio):',num2str(1-abs(sum(intT))/area2)]);
                    disp(['Error of integration of (spherical) Triangle (abs):',num2str(abs(sum(intT)-area2))]);
                end

            else % standard integration weights

                % standard integration weights on wedge
                [wint,~]  = this.W.ComputeIntegrationVector();

                % combine the integration weights
                int       = [wint,Tint];

                % compute the exact area of the whole shape
                hh   = abs(this.h);
                th   = pi - acos(hh/this.R);
                area = this.R^2*(th-0) + ...
                    this.R^2/2*(sin(2*0) - sin(2*th));
                this.Int  = int;
            end

            % display the error in the area if requested
            if(nargin == 2 && displayArea)

                if(area ~= 0)
                    disp(['BigSegment: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end

                disp(['BigSegment: Error of integration of area (absolute,area = 0): ',...
                    num2str(area-sum(this.Int))]);
            end

        end

        function Diff = ComputeDifferentiationMatrix(~)
            warning('Differentiation is not well-defined.')
            Diff = [];
        end

        %------------------------------------------------------------------
        % Plotting (generalized from Shape)
        %------------------------------------------------------------------

        function PlotGrid(this)
            % plot grid on wedge and triangle
            this.W.PlotGrid();
            hold on;
            this.T.PlotGrid();
        end

        function PlotGridLines(this)
            % plot grid lines on wedge and triangle
            this.W.PlotGridLines();
            hold on;
            this.T.PlotGridLines();
        end

        %------------------------------------------------------------------
        % Auxiliary functions
        %------------------------------------------------------------------

        function throwHRError(~)
            exc = MException('BigSegment:Constructor', 'Wall must intersect with disc');
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
            R = this.R;
            horizontal = this.horizontal; %#ok<*PROPLC>
            wallPos = this.Wall_Y;

            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * R + O1, 1.1 * R + O1];
            yLim = [-1.1 * R + O2, 1.1 * R + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot circle
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName', 'Circle and wall');

            % plot lines, and origin
            if(horizontal)
                plot(xLim, [wallPos, wallPos], '--k', 'HandleVisibility','off');
            else
                plot([wallPos, wallPos], yLim, '--k', 'HandleVisibility','off');
            end
            plot(O1,O2,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicator arrows
            % arrow for R
            p1 = [O1, O2];
            theta = -pi/2;
            p2 = [O1 + R * cos(theta), O2 + R * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R, p1(2) + dp(2)/2 - + 0.02 * R, 'R', ...
                'FontSize',20,'Color','b')

            % plot indicators for the two angles
            th1 = this.W.th1;
            th2 = this.W.th2;
            tth1 = linspace(0, mod(th1,2*pi), 100);
            tth2 = linspace(0, mod(th2,2*pi), 100);
            plot(R/4*cos(tth1) + O1, R/4*sin(tth1) + O2, ...
                '-m', 'DisplayName', 'th1, t2');
            plot(R/3*cos(tth2) + O1, R/3*sin(tth2) + O2, ...
                '-m', 'HandleVisibility','off');

            % line through origin to measure angles
            plot([O1,xLim(2)],[O2, O2], ':m', 'HandleVisibility','off');

            % plot indicators for the domain decomposition
            plot([O1, R*cos(th1) + O1],[O2, R*sin(th1) + O2],':k', ...
                'DisplayName','Domain split');
            plot([O1, R*cos(th2) + O1],[O2, R*sin(th2) + O2],':k', ...
                'HandleVisibility','off');

            % plot domain
            if(horizontal)
                tCircle = linspace(mod(th1,2*pi) - 2*pi, mod(th2,2*pi), 100);
            else
                tCircle = linspace(mod(th1,2*pi), mod(th2,2*pi), 100);
            end
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            plot([R*cos(th1) + O1, R*cos(th2) + O1],[R*sin(th1) + O2, R*sin(th2) + O2], ...
                '-r', 'HandleVisibility','off');

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

            geom.NW = [10;12];
            geom.NT = [10;12];
            geom.R = 1;
            geom.Origin = [-0.2,0.5];

            % geom.Wall_VertHor = 'horizontal';
            % geom.Wall_Y = geom.Origin(2) + 0.6;
            % geom.Wall = 'E';
            
            geom.Wall_VertHor = 'vertical';
            geom.Wall_Y = geom.Origin(1) + 0.6;
            geom.Wall = 'S';
            

        end

    end


end

