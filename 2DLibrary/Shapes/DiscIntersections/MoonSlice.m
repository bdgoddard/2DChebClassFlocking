classdef MoonSlice < M1SpectralSpectral
    % Part of a disc of radius R with origin (0,d) once the intersection
    % with a disc of radius R with origin (0,0) is remove and chopped off
    % along the horizontal lines through the intersection points.
    % Note: Origin is hard coded to 0.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radii of the two discs
    %  d - x-coordinate of the second disc

    properties
        R % radii of the two discs
        d % second disc is located at (0,d)
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        
        function this = MoonSlice(Geometry)

            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.R = Geometry.R;
            this.d = Geometry.d;
            this.Origin = [0;0]; % hard coded origin

            % ensure we're in polar coordinates
            this.polar = 'cart';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)
            
            % use the standard integration
            int = ComputeIntegrationVector@M1SpectralSpectral(this);

            % exact area
            % compute area of removed segments on left
            th   = 2*acos((this.d/2)/this.R);
            areaRemovedLeft = 2 * this.R^2/2*(th-sin(th));
            % compute area of removed segments on top and bottom
            xIntersect = this.d/2;
            yIntersect = sqrt(this.R^2 - xIntersect^2);
            th   = 2*acos(yIntersect/this.R);
            areaRemovedTB = 2 * this.R^2/2*(th-sin(th));

            % total area minus removed areas
            areaTotal = pi*this.R^2;
            area = areaTotal - areaRemovedLeft - areaRemovedTB;

            if(nargin == 2 && displayArea)
                if(area == 0)
                    disp(['MoonSlice: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['MoonSlice: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end
            end


        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------
        
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)

            % x and y coordinates of the (upper) intersection point
            xIntersect = this.d/2;
            yIntersect = sqrt(this.R^2 - xIntersect^2);

            % y2 map is between the two cut-off lines
            [y2_kv,dy2dx2] =  LinearMap(x2,-yIntersect,yIntersect);

            % determine left and right boundaries as a function of y2
            left = sqrt(this.R^2-y2_kv.^2);
            right = left + this.d;

            % linear map for each slice
            [y1_kv,dy1dx1] =  LinearMap(x1,left,right);

            % Jacobian
            if(nargout >= 3)
                J        = zeros(this.M,2,2);
                J(:,1,1) = dy1dx1;
                J(:,2,2) = dy2dx2;
            end

            if(nargout >= 4)
                dH1        = zeros(this.M,2,2);
            end

            if(nargout >= 4)
                dH2        = zeros(this.M,2,2);
            end
           
        end

        function [x1,x2] = CompSpace(this,y1,y2) %#ok
            % not yet implemented
            exc = MException('MoonSlice:CompSpace','not yet implemented');
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
            R = this.R; %#ok<*PROPLC>
            d = this.d;
            O1 = this.Origin(1); % should be hard coded to (0,0)
            O2 = this.Origin(2);
            OD1 = d;
            OD2 = 0;
            theta = acos((d/2)/this.R); % theta to the intersection point
            y2Min = O2 - R*sin(theta); % intersection point coordinates
            y2Max = O2 + R*sin(theta);
            y1Min = d/2;
            y1Max = OD1 + R*cos(theta);

            % set sensible axis limits
            xLim = [-1.1 * R + O1, 3.1 * R + O1];
            yLim = [-1.1 * R + O2, 1.1 * R + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circles, lines, and origins
            tCircle = linspace(0, 2*pi, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName', 'Circles and cutoff lines');
            plot(R*cos(tCircle) + OD1, R*sin(tCircle) + OD2, ...
                '--k', 'HandleVisibility','off');
            plot(xLim,[y2Min,y2Min], '--k', 'HandleVisibility','off');
            plot(xLim,[y2Max,y2Max], '--k', 'HandleVisibility','off');

            % lines

            % origins
            plot(O1,O2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','(0,0)')
            plot(OD1,OD2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','(d,0)')

            % plot indicators for R
            % R at (0,0)
            p1 = [O1, O2];
            thetaR = 7*pi/6;
            p2 = [O1 + R * cos(thetaR), O2 + R * sin(thetaR)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.03 * R, p1(2) + dp(2)/2 - 0.03 * R, 'R', ...
                'FontSize',20,'Color','b')

            % R at (d,0)
            p1 = [OD1, OD2];
            thetaR = 0;
            p2 = [OD1 + R * cos(thetaR), OD2 + R * sin(thetaR)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2), p1(2) + dp(2)/2 - 0.06 * R, 'R', ...
                'FontSize',20,'Color','b')

            % plot the boundary
            % arcs
            tCircle = linspace(-theta, theta, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            plot(R*cos(tCircle) + OD1, R*sin(tCircle) + OD2, ...
                '-r', 'HandleVisibility','off');
            % lines
            plot([y1Min,y1Max],[y2Min,y2Min], '-r', 'HandleVisibility','off');
            plot([y1Min,y1Max],[y2Max,y2Max], '-r', 'HandleVisibility','off');


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

            geom.N = [12;10];
            geom.R = 1;
            geom.d = 1.6;

        end

    end


end