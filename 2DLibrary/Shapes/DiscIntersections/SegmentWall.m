classdef SegmentWall < Segment
    % A segment of a disc of radius R, created by a
    % horizontal or vertical chord (determined by Wall_VertHor),
    % Wall_Y from the Origin.
    %
    % Inputs:
    %   N - [N1,N2] number of points in the two coordinates
    %   R - radius of the disc
    %   Wall_VertHor - 'vertical' or 'horizontal' orientation of wall
    %   Wall_Y - wall coordinate *not distance from wall*
    %   (Origin = [0,0]) - origin of the disc
    %   (sphere = false) - true implements 'spherical integration'
    %   (alpha) - fixed angle of rotation from the horizontal

    properties
        alpha = 0  % angle of rotation of points from horizontal
        horizontal % horizontal/vertical flag
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = SegmentWall(Geometry)
            
            % validate vertical/horizontal choice
            horizontal = strcmp(Geometry.Wall_VertHor,'horizontal');
            vertical   = strcmp(Geometry.Wall_VertHor,'vertical');
            if(~(horizontal || vertical))
                exc = MException('SegmentWall:Constructor', 'Wall_VertHor must be horizontal or vertical');
                throw(exc);
            end

            % determine the distance from the wall in the relevant
            % coordinate
            if(horizontal)
                Geometry.h = abs(Geometry.Wall_Y - Geometry.Origin(2));
            else % vertical
                Geometry.h = abs(Geometry.Wall_Y - Geometry.Origin(1));
            end

            % use the Segment constructor to create a segment in the upper
            % half plane distance |h| from the wall
            this@Segment(Geometry);

            % determine if the wall is horizontal
            this.horizontal = horizontal;

            % determine which half plane we should actually be in
            if(horizontal && (Geometry.Wall_Y > Geometry.Origin(2)))
                % upper half plane
                this.alpha = 0;
            elseif(horizontal && (Geometry.Wall_Y <= Geometry.Origin(2)))
                % lower half plane
                this.alpha = -pi;
            elseif(vertical && (Geometry.Wall_Y > Geometry.Origin(1)))
                % right half plane
                this.alpha = -pi/2;
            elseif(vertical && (Geometry.Wall_Y <= Geometry.Origin(1)))
                % left half plane
                this.alpha = pi/2;
            end

            % can override to produce segment at an arbitrary angle
            if(isfield(Geometry,'alpha'))
                this.alpha = Geometry.alpha;
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % map is not invertible at the corners
            this.invertibleMap = false; 

        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J] = PhysSpace(this,x1,x2)
            % use the default points from Segment
            [y1_kv,y2_kv,J]  = PhysSpace@Segment(this,x1,x2);
            % rotate to the correct angle
            [y1_kv,y2_kv,JR] = RotateAroundOrigin(y1_kv,y2_kv,this.alpha);
            % compute the resulting Jacobian
            J                = MultMatr(JR,J);
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
            % determine signed h
            if(horizontal)
                h = this.h * cos(this.alpha);
            else
                h = -this.h * sin(this.alpha);
            end
            O1 = this.Origin(1);
            O2 = this.Origin(2);
            theta = asin(h/R);

            % set sensible axis limits
            xLim = [-1.1 * R + O1, 1.1 * R + O1];
            yLim = [-1.1 * R + O2, 1.1 * R + O2];

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
            if(horizontal)
                plot(xLim,[h + O2, h + O2], ...
                    '--k',  'HandleVisibility','off');
                % line through origin to measure h from
                plot(xLim,[O2, O2], ':k', 'HandleVisibility','off');
            else
                plot([h + O1, h + O1], yLim, ...
                    '--k',  'HandleVisibility','off');
                % line through origin to measure h from
                plot([O1, O1], yLim, ':k', 'HandleVisibility','off');
            end
            plot(O1,O2,'ok','MarkerSize',15,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicators
            % indicator for h
            if(horizontal)
                p1 = [O1 - R/2, O2];
                p2 = [O1 - R/2, O2 + h];
                dp = p2 - p1;
                text((p1(1) + dp(1)/2) + 0.02 * R, p1(2) + dp(2)/2, 'h', ...
                    'FontSize',20,'Color','b')
            else
                p1 = [O1, O2 - R/2];
                p2 = [O1 + h, O2 - R/2];
                dp = p2 - p1;
                text((p1(1) + dp(1)/2), p1(2) + dp(2)/2 + 0.05 * R, 'h', ...
                    'FontSize',20,'Color','b')
            end
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');

            % indicator for R
            p1 = [O1, O2];
            thetaR = pi/6;
            p2 = [O1 + R * cos(thetaR), O2 + R * sin(thetaR)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R, p1(2) + dp(2)/2 - + 0.02 * R, 'R', ...
                'FontSize',20,'Color','b')


            % plot the boundary
            % circle
            if(horizontal)
                tCircle = linspace(theta, pi-theta, 100);
            else
                tCircle = linspace(theta - pi/2, pi/2-theta, 100);
            end
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            % line
            if(horizontal)
                y1Min = O1 - R*cos(theta);
                y1Max = O1 + R*cos(theta);
                plot([y1Min,y1Max],[O2+h,O2+h],'r','HandleVisibility','off');
            else
                y2Min = O2 - R*cos(theta);
                y2Max = O2 + R*cos(theta);
                plot([O1+h,O1+h],[y2Min,y2Max],'r','HandleVisibility','off');
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
            % geom.Wall_VertHor = 'vertical';
            geom.Wall_VertHor = 'horizontal';
            geom.Wall_Y = 0.2;
            geom.R = 1;
            geom.Origin = [0,0];
            
        end

    end


end