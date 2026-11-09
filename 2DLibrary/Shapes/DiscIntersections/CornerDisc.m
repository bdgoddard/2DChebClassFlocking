classdef CornerDisc < handle
    % Part a disc of radius R, centred at Origin, created by it being cut
    % off by the corner of an abritrarily large square, with the corner at
    % CornerPos.
    %
    % Note that the corner must lie inside the disc.  Otherwise try
    % DoubleCutDisc.
    %
    % This is formed from a triangle and segment.
    %
    % Inputs:
    %   NT - [NT1,NT2] number of points for the triangle in the two coordinates
    %   NS - [NS1,NS2] number of points for the segment in the two
    %   coordinates
    %   R - radius of the disc
    %   (Corner = 'SW') - one of 'NE', 'SE', 'SW', 'NW' denoting the compass direction of
    %   the corner from the centre of the square
    %   (CornerPos = [0,0]) - [CP1,CP2] the position of the corner
    %   (Origin = [0,0]) - origin of the disc

    properties
        R  % radius of disc
        Origin = [0;0]; % origin of disc
        CornerPos = [0;0]; % position of corner
        Corner = 'SW'; % orientation of corner
        S % segment
        T % triangle
        Int % integration vector
        Pts % points structure
    end

    methods

        function this = CornerDisc(Geometry)

            % set parameters
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            if(isfield(Geometry,'CornerPos'))
                this.CornerPos = Geometry.CornerPos;
            end

            if(isfield(Geometry,'Corner'))
                this.Corner = Geometry.Corner;
            end

            % compute geometry for segment:
            this.R      = Geometry.R;

            % offset of Corner and Origin
            dOC1 = this.Origin(1) - this.CornerPos(1);
            dOC2 = this.Origin(2) - this.CornerPos(2);

            % check that the corner lies inside the disc
            if(this.R^2 <= (dOC1^2 + dOC2^2))
                exc = MException('CornerDisc:Constructor', ...
                    'Corner must lie inside Disc.  Try using DoubelCutDisc.');
                throw(exc);
            end

            % points (dw1,CornerPos(2)), (CornerPos(1),dw2)
            % at which the circle intersects with corner relative to the origin
            dw1 = sqrt(this.R^2 - dOC2^2);
            dw2 = sqrt(this.R^2 - dOC1^2);

            % points at which the circle intersects the corner in physical
            % space
            % and the angle at which the Segment sits
            switch(this.Corner)
                case 'SW'
                    yi1 = this.Origin(1) + dw1;
                    yi2 = this.Origin(2) + dw2;
                    theta = -atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
                case 'NW'
                    yi1 = this.Origin(1) + dw1;
                    yi2 = this.Origin(2) - dw2;
                    theta = pi-atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
                case 'NE'
                    yi1 = this.Origin(1) - dw1;
                    yi2 = this.Origin(2) - dw2;
                    theta = pi -atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
                case 'SE'
                    yi1 = this.Origin(1) - dw1;
                    yi2 = this.Origin(2) + dw2;
                    theta = -atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
            end

            % compute points for the Segment
            segSh.N = Geometry.NS;
            segSh.R = this.R;
            % midpoint of the line y1p-y2p
            y1p = 1/2*(yi1 + this.CornerPos(1));
            y2p = 1/2*(yi2 + this.CornerPos(2));
            % height of the segment above the imaginary wall
            segSh.h = sqrt( (this.Origin(1) - y1p)^2 + (this.Origin(2) - y2p)^2 );
            this.S  = Segment(segSh);
            % rotate the result by theta and shift the origin to [y1p,y2p]
            ShiftRotateSegment(this,segSh.h,theta,[y1p;y2p]);

            % compute points for the Triangle:
            % define the corners
            triSh.Y = [this.CornerPos(1), this.CornerPos(2);
                this.CornerPos(1), yi2;
                yi1, this.CornerPos(2)];
            triSh.N = Geometry.NT;
            this.T = Triangle(triSh);

            % concatenate the points
            this.Pts  = struct('y1_kv',[this.S.Pts.y1_kv; this.T.Pts.y1_kv],...
                'y2_kv',[this.S.Pts.y2_kv; this.T.Pts.y2_kv]);

        end

        %------------------------------------------------------------------
        % Shape-like functions
        %------------------------------------------------------------------

        % Cartesian points
        function ptsCart = GetCartPts(this)
            ptsCart = this.Pts;
        end

        function [int,area] = ComputeIntegrationVector(this)

            % integration weights for the triangle and segment
            [Tint,areaT]  = this.T.ComputeIntegrationVector();
            [Sint,areaS]  = this.S.ComputeIntegrationVector();

            % concatenate the integrtion vectors
            int  = [Sint,Tint];

            % explicit area
            area = areaT + areaS;

            % assign to Int
            this.Int = int;

        end

        function Diff = ComputeDifferentiationMatrix(~)
            warning('Differentiation is not well-defined.')
            Diff = [];
        end

        %------------------------------------------------------------------
        % Auxiliary functions
        %------------------------------------------------------------------

        function ShiftRotateSegment(this,h,th,yp)
            % shift and rotate the segment to the correct position

            y1_kv = this.S.Pts.y1_kv;
            y2_kv = this.S.Pts.y2_kv - h;

            y1_kvs = cos(th)*y1_kv - sin(th)*y2_kv + yp(1);
            y2_kvs = sin(th)*y1_kv + cos(th)*y2_kv + yp(2);

            this.S.Pts.y1_kv = y1_kvs;
            this.S.Pts.y2_kv = y2_kvs;

        end

        %------------------------------------------------------------------
        % Plotting (generalized from Shape)
        %------------------------------------------------------------------

        function PlotGrid(this)
            % plot grid on wedge and triangle
            this.S.PlotGrid();
            hold on;
            this.T.PlotGrid();
        end

        function PlotGridLines(this)
            % plot grid lines on wedge and triangle
            this.S.PlotGridLines();
            hold on;
            this.T.PlotGridLines();
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
            R = this.R; %#ok<*PROPLC>
            yi1 = this.T.Y(3,1);
            yi2 = this.T.Y(2,2);

            O1 = this.Origin(1);
            O2 = this.Origin(2);
            CP1 = this.CornerPos(1);
            CP2 = this.CornerPos(2);

            % set sensible axis limits
            xLim = [-1.1 * R + O1, 1.1 * R + O1];
            yLim = [-1.1 * R + O2, 1.1 * R + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circle
            tCircle = linspace(0, 2*pi, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName', 'Circle and corner');

            % plot the corner
            switch(this.Corner)
                case 'SW'
                    plot([CP1,CP1],[CP2,yLim(2)],'--k','HandleVisibility','off');
                    plot([CP1,xLim(2)],[CP2,CP2],'--k','HandleVisibility','off');
                case 'NW'
                    plot([CP1,CP1],[yLim(1),CP2],'--k','HandleVisibility','off');
                    plot([CP1,xLim(2)],[CP2,CP2],'--k','HandleVisibility','off');
                case 'NE'
                    plot([CP1,CP1],[yLim(1),CP2],'--k','HandleVisibility','off');
                    plot([xLim(1),CP1],[CP2,CP2],'--k','HandleVisibility','off');
                case 'SE'
                    plot([CP1,CP1],[CP2,yLim(2)],'--k','HandleVisibility','off');
                    plot([xLim(1),CP1],[CP2,CP2],'--k','HandleVisibility','off');
            end

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',15,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicator arrows
            % arrow for R
            p1 = [O1, O2];
            theta = pi/6;
            p2 = [O1 + R * cos(theta), O2 + R * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R, p1(2) + dp(2)/2 - + 0.02 * R, 'R', ...
                'FontSize',20,'Color','b')

            % plot indicators for the domain decomposition
            plot([CP1, yi1],[yi2, CP2],':k', ...
                'DisplayName','Domain split');

            % plot boundary
            switch(this.Corner)
                case 'SW'
                    plot([CP1,CP1],[CP2,yi2],'-r','DisplayName', 'Boundary');
                    plot([CP1,yi1],[CP2,CP2],'-r','HandleVisibility','off');
                    th1 = pi/2 - asin((CP1-O1)/R);
                    th2 = asin((CP2-O2)/R);
                case 'NW'
                    plot([CP1,CP1],[yi2,CP2],'-r','DisplayName', 'Boundary');
                    plot([CP1,yi1],[CP2,CP2],'-r','HandleVisibility','off');
                    th1 = -pi/2 + asin((CP1-O1)/R);
                    th2 = asin((CP2-O2)/R);
                case 'NE'
                    plot([CP1,CP1],[yi2,CP2],'-r','DisplayName', 'Boundary');
                    plot([yi1,CP1],[CP2,CP2],'-r','HandleVisibility','off');
                    th1 = pi/2 + acos((CP2-O2)/R);
                    th2 = 3*pi/2 + asin((CP1-O1)/R);
                case 'SE'
                    plot([CP1,CP1],[CP2,yi2],'-r','DisplayName', 'Boundary');
                    plot([yi1,CP1],[CP2,CP2],'-r','HandleVisibility','off');
                    th1 = pi/2 - asin((CP1-O1)/R);
                    th2 = pi/2 + acos((CP2-O2)/R);                    
            end
            % plot the arc
            tCircle = linspace(th1, th2, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
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

            geom.NS = [10;12];
            geom.NT = [10;12];
            geom.R = 1;
            geom.CornerPos = [-0.4,0.3];
            geom.Origin = [-0.1,-0.2];
            geom.Corner = 'NW';

        end

    end


end

