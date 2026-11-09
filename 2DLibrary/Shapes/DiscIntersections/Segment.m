classdef Segment < M1SpectralSpectral
    % A segment of a disc of radius R, created by a
    % horizontal chord at (positive) distance h above the origin.
    % To be used when h < R; otherwise use BigSegment.
    %
    % Inputs:
    %   N - [N1,N2] number of points in the two coordinates
    %   R - radius of the disc
    %   h - height of chord above the origin
    %   (Origin = [0,0]) - origin of the disc
    %   (sphere = false) - true implements 'spherical integration'

    properties
        h  % height of chord above the Origin
        R  % radius of the circle
        sphere = false  % integrate as if it's a cap of a sphere
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = Segment(Geometry)

            % validate h and R
            if(Geometry.h > Geometry.R)
                exc = MException('Segment:Constructor', 'Requires h <= R; use BigSegment');
                throw(exc);
            end

            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.R = Geometry.R;
            this.h = abs(Geometry.h); % note we enforce h>=0

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % set geometry flags
            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;
            end

            % we're in Cartesian coordinates
            this.polar = 'cart';

            % map is not invertible at the corners
            this.invertibleMap = false; 

        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)

            int = ComputeIntegrationVector@M1SpectralSpectral(this);

            % compute exact area
            if(this.sphere)

                % rescale int for spherically-weighted integration
                y1s  = this.Pts.y1_kv;
                y2s  = this.Pts.y2_kv;
                int  = 2*int.*real(sqrt(this.R^2 - y1s.^2 - y2s.^2))';
                this.Int = int; % override integration vector

                % explicit weighted area of a segment
                ht   = this.R - this.h;
                area = pi*ht^2/3*(3*this.R - ht);

            else % standard integration

                % explicit area of a segment
                th   = 2 * acos( this.h / this.R );
                area = this.R^2 / 2 * ( th-sin(th) );

            end

            % display the error in the area if requested
            if(nargin == 2 && displayArea)

                if(area == 0)
                    disp(['Segment: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['Segment: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end

            end

        end % ComputeIntegrationVector

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)

            % Note that x1 and x2 are x1_kv and x2_kv

            % radius of circle and chord height above the Origin
            R  = this.R; %#ok<*PROPLC>
            h  = this.h;

            n  = length(x1);  % number of points used to initialise matrices

            % (positive) y1 value at intersection of chord and circle
            C1      = sqrt(R^2-h^2);
            y1_kv   = C1 * x1;  % scale [-1,1] to the base of the Segment

            % compute the height of the segment at each y1
            % height of circle is sqrt(R^2-y1_kv.^2)
            % and subtract h
            g       = sqrt(R^2-y1_kv.^2) - h;

            % derivatives of g (explicit calculation)
            dgdx1   = -C1^2 * x1 ./ (sqrt(R^2-y1_kv.^2));
            ddgddx1 = -C1^2 * R^2 ./ ((R^2-y1_kv.^2).^(3/2));

            % scale [-1,1] to the height of the segment (g)
            % and shift upwards by h
            y2_kv   = ((1+x2)/2).*g + h;

            % Jacobian (explicit calculation)
            % Note that differentiation is not well defined as multiple
            % points get collapsed onto the corners
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = C1 * ones(n,1);
                J(:,2,1) = ((1+x2)/2) .* dgdx1;
                J(:,2,2) = g/2;
            end

            % Hessians (explicit calculation)
            if(nargout >= 4)
                dH1        = zeros(n,2,2);
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);
                dH2(:,1,1) = ddgddx1 .* (1+x2)/2;
                dH2(:,1,2) = dgdx1/2;
                dH2(:,2,1) = dgdx1/2;
            end

        end

        function [x1,x2] = CompSpace(this,y1,y2) %#ok
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
            R = this.R;
            h = this.h;
            O1 = this.Origin(1);
            O2 = this.Origin(2);
            theta = asin(h/R);

            % set sensible axis limits
            xLim = [-1.4 * R + O1, 1.4 * R + O1];
            yLim = [-1.4 * R + O2, 1.4 * R + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;

            % plot the circle, lines, and origin
            tCircle = linspace(0, 2*pi, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName', 'Circle and wall');
            plot(xLim,[h + O2, h + O2], ...
                '--k', 'HandleVisibility','off');
            plot(O1,O2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','Origin')

            % line through origin to measure h from
            plot(xLim,[O2, O2],':k', 'HandleVisibility','off');

            % plot indicators
            % indicator for h
            p1 = [O1 - R/2, O2];
            p2 = [O1 - R/2, O2 + h];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R, p1(2) + dp(2)/2, 'h', ...
                'FontSize',20,'Color','b')
            % indicator for R
            p1 = [O1, O2];
            thetaR = -pi/6;
            p2 = [O1 + R * cos(thetaR), O2 + R * sin(thetaR)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.04 * R, p1(2) + dp(2)/2 + 0.04 * R, 'R', ...
                'FontSize',20,'Color','b')

            % plot the boundary
            % circle
            tCircle = linspace(theta, pi-theta, 100);
            hold on
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            % line
            y1Min = O1 - R*cos(theta);
            y1Max = O1 + R*cos(theta);
            plot([y1Min,y1Max],[O2+h,O2+h],'r','HandleVisibility','off');


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

            geom.N = [10;12];
            geom.h = 0.3;
            geom.R = 1;
            geom.Origin = [-0.5,0.2];

        end

    end


end % classdef