classdef DoubleCutDisc < handle
    % Part a disc of radius R, centred at Origin, created by it being cut
    % off by the corner of an abritrarily large square, with the corner at
    % CornerPos.
    %
    % Note that the corner must lie outside the disc.  Otherwise try
    % CornerDisc.
    %
    % This is formed from two triangles and two wedges
    %
    % Inputs:
    %   NT - [NT1,NT2] number of points for the triangles in the two coordinates
    %   NW - [NS1,NS2] number of points for the wedges in the two coordinates
    %   R - radius of the disc
    %   (Corner = 'SW') - one of 'NE', 'SE', 'SW', 'NW' denoting the compass direction of
    %   the corner from the centre of the square
    %   (CornerPos = [0,0]) - [CP1,CP2] the position of the corner
    %   (Origin = [0,0]) - origin of the disc


    properties
        R % radius of disc
        Origin = [0;0]; % origin of disc
        CornerPos = [0;0]; % position of corner
        Corner    = 'SW'; % orientation of corner
        T1 % triangle 1
        T2 % triangle 1
        W1 % wedge 1
        W2 % wedge 2
        Int % integration vector
        Pts % points structure
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = DoubleCutDisc(Geometry)

            if(isfield(Geometry,'Origin'))
                this.Origin     = Geometry.Origin;
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

            % check that the corner lies outside the disc
            if(this.R^2 > (dOC1^2 + dOC2^2))
                exc = MException('DoubelCutDisc:Constructor', ...
                    'Corner must lie outside Disc.  Try using DoubelCutDisc.');
                throw(exc);
            end

            % points (dw1,CornerPos(2)), (CornerPos(1),dw2)
            % at which the circle intersects with corner relative to the origin
            dw1 = sqrt(this.R^2 - dOC2^2);
            dw2 = sqrt(this.R^2 - dOC1^2);

            % angles describing the domain decomposition
            theta1 = atan( dw2 / dOC1 );
            theta2 = atan( dw1 / dOC2 );

            % make sure that the angles are correct based on the
            % orientation of the corner
            switch(this.Corner)
                case 'SW'
                    wedgeSh1.th1 = pi + theta1;
                    wedgeSh1.th2 = 3*pi/2 - theta2;
                    wedgeSh2.th1 = -pi/2 + theta2;
                    wedgeSh2.th2 =  pi - theta1;
                case 'NW'
                    theta2 = -theta2;
                    wedgeSh1.th1 = pi/2 + theta2;
                    wedgeSh1.th2 = pi - theta1;
                    wedgeSh2.th1 = -pi + theta1;
                    wedgeSh2.th2 = pi/2 - theta2;
                case 'NE'
                    theta1 = -theta1;
                    theta2 = -theta2;
                    wedgeSh1.th1 = theta1;
                    wedgeSh1.th2 = pi/2 - theta2;
                    wedgeSh2.th1 = pi/2 + theta2;
                    wedgeSh2.th2 = 2*pi-theta1;
                case 'SE'
                    theta1 = -theta1;
                    wedgeSh1.th1 = 3*pi/2 + theta2;
                    wedgeSh1.th2 = 2*pi - theta1;
                    wedgeSh2.th1 = theta1;
                    wedgeSh2.th2 = 3*pi/2 - theta2;
            end

            % compute points for triangles
            triSh1.Y = [this.Origin(1), this.Origin(2);
                this.CornerPos(1),  dw2 + this.Origin(2);
                this.CornerPos(1), -dw2 + this.Origin(2)];
            triSh1.N = Geometry.NT;
            this.T1 = Triangle(triSh1);

            triSh2.Y = [this.Origin(1), this.Origin(2);
                dw1 + this.Origin(1), this.CornerPos(2);
                -dw1 + this.Origin(1), this.CornerPos(2)];
            triSh2.N = Geometry.NT;
            this.T2 = Triangle(triSh2);

            % compute points for the wedges
            wedgeSh1.R  = this.R;
            wedgeSh1.N  = Geometry.NW;
            wedgeSh1.Origin = this.Origin;
            this.W1 = Wedge(wedgeSh1);
            % this.W1.Origin = this.Origin; % set the correct origin

            wedgeSh2.R  = this.R;
            wedgeSh2.N  = Geometry.NW;
            wedgeSh2.Origin = this.Origin;
            this.W2 = Wedge(wedgeSh2);

            % concatenate the points
            this.Pts  = struct('y1_kv',[this.T1.GetCartPts.y1_kv; ...
                this.T2.GetCartPts.y1_kv; ...
                this.W1.GetCartPts.y1_kv; ...
                this.W2.GetCartPts.y1_kv],...
                'y2_kv',[this.T1.GetCartPts.y2_kv; ...
                this.T2.GetCartPts.y2_kv; ...
                this.W1.GetCartPts.y2_kv; ...
                this.W2.GetCartPts.y2_kv]);

        end


        %------------------------------------------------------------------
        % Shape-like functions
        %------------------------------------------------------------------

        % Cartesian points
        function ptsCart = GetCartPts(this)
            ptsCart = this.Pts;
        end

        function [int,area] = ComputeIntegrationVector(this)

            % integration weights for the triangles and wedges
            [T1int,area1]  = this.T1.ComputeIntegrationVector();
            [T2int,area2]  = this.T2.ComputeIntegrationVector();
            [W1int,area3]  = this.W1.ComputeIntegrationVector();
            [W2int,area4]  = this.W2.ComputeIntegrationVector();

            % concatenate the integrtion vectors
            int  = [T1int,T2int,W1int,W2int];

            % explicit area
            area = area1 + area2 + area3 + area4;

            % assign to Int
            this.Int = int;

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
            this.T1.PlotGrid();
            hold on;
            this.T2.PlotGrid();
            this.W1.PlotGrid();
            this.W2.PlotGrid();
        end

        function PlotGridLines(this)
            % plot grid lines on wedge and triangle
            this.T1.PlotGridLines();
            hold on;
            this.T2.PlotGridLines();
            this.W1.PlotGridLines();
            this.W2.PlotGridLines();
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
            plot(O1,O2,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','Origin')

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
            T1Y = this.T1.Y;
            T2Y = this.T2.Y;
            plot([O1,T1Y(2,1)],[O2,T1Y(2,2)],':k','DisplayName','Domain split');
            plot([O1,T1Y(3,1)],[O2,T1Y(3,2)],':k','HandleVisibility','off');
            plot([O1,T2Y(2,1)],[O2,T2Y(2,2)],':k','HandleVisibility','off');
            plot([O1,T2Y(3,1)],[O2,T2Y(3,2)],':k','HandleVisibility','off');

            % plot the boundary
            % edges
            plot([CP1,CP1],[T1Y(2,2),T1Y(3,2)],'-r','DisplayName','Boundary');
            plot([T2Y(2,1),T2Y(3,1)],[CP2,CP2],'-r','HandleVisibility','off');
            % arcs
            switch(this.Corner)
                case 'SW'
                    th1 = acos((CP1 - O1)/R);
                    th2 = - acos((CP1 - O1)/R);
                    th3 = -pi/2 - acos(-(CP2 - O2)/R);
                    th4 = -pi/2 + acos(-(CP2 - O2)/R);
                case 'NW'
                    th1 = acos((CP1 - O1)/R);
                    th2 = - acos((CP1 - O1)/R);
                    th3 = 3*pi/2 -  acos(-(CP2 - O2)/R);
                    th4 = -pi/2 + acos(-(CP2 - O2)/R);
                case 'NE'
                    th1 = 2*pi - acos((CP1 - O1)/R);
                    th2 = acos((CP1 - O1)/R);
                    th3 = 3*pi/2 -  acos(-(CP2 - O2)/R);
                    th4 = -pi/2 + acos(-(CP2 - O2)/R);
                case 'SE'
                    th1 = acos((CP1 - O1)/R);
                    th2 = - acos((CP1 - O1)/R);
                    th3 = 3*pi/2 -  acos(-(CP2 - O2)/R);
                    th4 = -pi/2 + acos(-(CP2 - O2)/R);                    
            end
            tCircle = linspace(th2, th4, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'HandleVisibility','off');
            tCircle = linspace(th1, th3, 100);
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

            geom.NW = [10;12];
            geom.NT = [10;12];
            geom.R = 1;
            geom.Origin = [0.7,0.85];

            % geom.CornerPos = [geom.Origin(1) - 0.9*geom.R, geom.Origin(2) - 0.75*geom.R]; % SW
            % geom.Corner = 'SW';

            % geom.CornerPos = [geom.Origin(1) - 0.9*geom.R, geom.Origin(2) + 0.75*geom.R]; % NW
            % geom.Corner = 'NW';

            % geom.CornerPos = [geom.Origin(1) + 0.9*geom.R, geom.Origin(2) + 0.75*geom.R]; % NE
            % geom.Corner = 'NE';

            geom.CornerPos = [geom.Origin(1) + 0.9*geom.R, geom.Origin(2) - 0.75*geom.R]; % NE
            geom.Corner = 'SE';

        end

    end

end

