classdef Moon < Polar_M1SpectralSpectral
    % Part of a disc of radius R with origin (0,d) once the intersection
    % with a disc of radius R with origin (0,0) is removed.
    % Note: must have d < R.
    % Note: Origin is hard coded to 0.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radii of the two discs
    %  d - x-coordinate of the second disc

    properties
        R  % radii of the two discs
        d  % second disc is located at (0,d)
        theta  % construction angle
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        
        function this = Moon(Geometry)

            % validate d and R
            if(Geometry.d >= Geometry.R)
                exc = MException('Moon:Constructor', 'Requires d < R.');
                throw(exc);
            end

            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.R    = Geometry.R;
            this.d    = Geometry.d;
            this.Origin = [0;0]; % force the Origin to be [0;0];

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % map is not invertible at the corners
            this.invertibleMap = false; 

        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)
            int = ComputeIntegrationVector@Polar_M1SpectralSpectral(this);

            % compute area of removed segments
            th   = 2*acos((this.d/2)/this.R);
            areaRemoved = 2 * this.R^2/2*(th-sin(th));
            
            % compute total area
            areaTotal = pi*this.R^2;
            area = areaTotal - areaRemoved;

            % display the error in the area if requested
            if(nargin == 2 && displayArea)
                if(area == 0)
                    disp(['Moon: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['Moon: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end
            end
        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------
        
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            % inputs are already kron products

            % linear in theta
            xIntersect = this.d/2; % discs are are at (0,0) and (0,d)
            this.theta = acos(xIntersect/this.R); % theta to the intersection point
            [y2_kv,J2] = LinearMap(x2,-this.theta,this.theta); % angle at origin

            % 'Inner' edge
            RIn = this.R*ones(size(x1));

            % compute ROut by solving a quadratic (by hand) given that
            % 1) The point (x,y) lies on the circle of radius R centred
            % around (d,0)
            % 2) tan(theta) = y/x

            a = 1 + tan(y2_kv).^2;
            b = -2*this.d;
            c = this.d^2 - this.R^2;
            xOut = ( - b + sqrt(b.^2 - 4 * a * c) ) ./ (2 * a);
            yOut = xOut.*tan(y2_kv);
            ROut = sqrt(xOut.^2 + yOut.^2);

            % standard Chebyshev map from RIn to ROut; different for
            % different theta = y2
            [y1_kv,J1] = LinearMap(x1,RIn,ROut);

            if(nargout >= 3)
                J        = zeros(this.M,2,2);
                J(:,1,1) = J1; %(ROut - RIn)/2;
                J(:,2,2) = J2; %(pi-2*theta)/2;
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
            exc = MException('Moon:CompSpace','not yet implemented');
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
            theta = this.theta;
            O1 = this.Origin(1); % should be hard coded to (0,0)
            O2 = this.Origin(2);
            OD1 = d;
            OD2 = 0;

            % set sensible axis limits
            xLim = [-1.1 * R + O1, 2.1 * R + O1];
            yLim = [-1.1 * R + O2, 1.1 * R + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circles, and origins
            tCircle = linspace(0, 2*pi, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName', 'Circle of radius R at (0,0)');
            plot(R*cos(tCircle) + OD1, R*sin(tCircle) + OD2, ...
                '--k', 'DisplayName', 'Circle of radius R at (d,0)');

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
            text((p1(1) + dp(1)/2), p1(2) + dp(2)/2 - 0.05 * R, 'R', ...
                'FontSize',20,'Color','b')

            % plot the boundary
            tCircle = linspace(-theta, theta, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            tCircle = linspace(pi - theta, -pi + theta, 100);
            plot(R*cos(tCircle) + OD1, R*sin(tCircle) + OD2, ...
                '-r', 'HandleVisibility','off');

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
            geom.d = 0.8;

        end

    end
    
end