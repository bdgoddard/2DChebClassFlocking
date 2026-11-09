classdef SlitMinusDisk < M1SpectralSpectral
    % Part of a half-infinite slit (either 'Left' or 'Right') cut off at 
    % one end by a disc of radius R centred at OriginDisk.  The finite end 
    % of the strip is at OriginDisk(1).
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radius of the disc to be removed
    %  L1 - mapping parameter for radial coordinate; half the points lie in
    %       [R,R+L1]
    %  y2Min, y2Max - positions of the upper and lower walls
    %  (LeftRight = 'Right') - 'Left' or 'Right' for direction of slit
    %  (OriginDisk = [0,0]) - origin of the disc

    properties
        R  % radius of disc
        OriginDisk = [0,0] % origin of disc
        y2Min % bottom wall position
        y2Max % top wall position
        L1 % mapping parameter - (roughly) hald the points are in 
           % [OriginDisk(1),OriginDisk(1)+L1]
        LeftRight = 'Right'; % 'Left' or 'Right' orientation of slit
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = SlitMinusDisk(Geometry)

            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.OriginDisk = Geometry.OriginDisk;
            this.y2Min      = Geometry.y2Min;
            this.y2Max      = Geometry.y2Max;
            this.LeftRight  = Geometry.LeftRight;

            % ensure we're in Cartesian coordinates
            this.polar = 'cart';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
            
        end


        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@M1SpectralSpectral(this);
            area = inf;
        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------
        
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)

            % quotient map from 0 to inf
            [y1,dy1]    = QuotientMap(x1,this.L1,0,inf);
            % linear map across the strip
            [y2_kv,dy2] = LinearMap(x2,this.y2Min,this.y2Max);

            % find the left hand edge for each y2
            y1Shift    = sqrt(this.R^2 - (y2_kv - this.OriginDisk(2)).^2);

            % shift/flip the points depending on the orientation
            if(strcmp(this.LeftRight,'Right'))
                y1_kv  = y1 + this.OriginDisk(1) + real(y1Shift);
            elseif(strcmp(this.LeftRight,'Left'))
                y1_kv  = -y1 + this.OriginDisk(1) - real(y1Shift);
                dy1    = - dy1;
            else
                error('SlitMinusDisk:choose "Left" or "Right" for this.LeftRight.');
            end

            % Jacobian and Hessian
            n = this.N1*this.N2;
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = dy1;
                %J(:,2,1) = ((1+x2)/2).*dgdx1;
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
            y2Min = this.y2Min;
            y2Max = this.y2Max;
            L2 = y2Max - y2Min;
            OD1 = this.OriginDisk(1);
            OD2 = this.OriginDisk(2);
            right = strcmp(this.LeftRight,'Right');
            insideTop = ( abs( y2Max - OD2 ) > R );
            insideBottom = ( abs( y2Min - OD2 ) > R );

            % set sensible axis limits
            if(right)
                xLim = [-1.1 * L1 +  OD1, 1.1 * 3*L1 + OD1];
            else
                xLim = [-1.1 * 3*L1 + OD1, 1.1 * L1 + OD1];
            end
            % deal with the case when the strip is narrow
            y2LimMin = min(y2Min - 0.5*L2,OD2 - R - 0.5*L2);
            y2LimMax = max(y2Max + 0.5*L2,OD2 + R + 0.5*L2);
            yLim = [y2LimMin,y2LimMax];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator line for L1
            if(right)
                plot([OD1+L1,OD1+L1],[y2Min,y2Max],'--m','DisplayName','OriginDisk(1)+L1');
            else
                plot([OD1-L1,OD1-L1],[y2Min,y2Max],'--m','DisplayName','OriginDisk(1)-L1');
            end

            % plot circle
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + OD1, R*sin(tCircle) + OD2, ...
                '--k', 'DisplayName','Disc and slit');

            % plot walls
            if(right)
                plot([OD1,xLim(2)],[y2Max,y2Max],'--k', 'HandleVisibility','off');
                plot([OD1,xLim(2)],[y2Min,y2Min],'--k', 'HandleVisibility','off');
            else
                plot([xLim(1),OD1],[y2Max,y2Max],'--k', 'HandleVisibility','off');
                plot([xLim(1),OD1],[y2Min,y2Min],'--k', 'HandleVisibility','off');
            end
            plot([OD1,OD1],[y2Min,y2Max],'--k', 'HandleVisibility','off');

            % plot the origin
            plot(OD1,OD2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','OriginDisk')

            % plot boundary
            % horizontal lines
            % deal with whether the disc is inside the slit
            if(~insideBottom)
                if(right)
                    xLimB = OD1 + sqrt( R^2 - (OD2 - y2Min)^2 );
                else
                    xLimB = OD1 - sqrt( R^2 - (OD2 - y2Min)^2 );
                end
            else
                xLimB = OD1;
            end
            if(~insideTop)
                if(right)
                    xLimT = OD1 + sqrt( R^2 - (OD2 - y2Max)^2 );
                else
                    xLimT = OD1 - sqrt( R^2 - (OD2 - y2Max)^2 );
                end
            else
                xLimT = OD1;
            end
            % deal with the orientation
            if(right)
                xRangeB = [xLimB,xLim(2)];
                xRangeT = [xLimT,xLim(2)];
            else
                xRangeB = [xLim(1),xLimB];
                xRangeT = [xLim(1),xLimT];
            end
            % plot the lines
            plot(xRangeB,[y2Min,y2Min],'-r', 'DisplayName', 'Boundary');
            plot(xRangeT,[y2Max,y2Max],'-r', 'HandleVisibility','off');

            % vertical lines)
            if(insideBottom)
                plot([OD1,OD1],[y2Min,OD2-R],'-r', 'HandleVisibility','off');
            end
            if(insideTop)
                plot([OD1,OD1],[OD2+R,y2Max],'-r', 'HandleVisibility','off');
            end
            
            % curved edge
            if(~insideTop)
                if(right)
                    thetaT = -asin((OD2 - y2Max)/R);
                else
                    thetaT = pi + asin((OD2 - y2Max)/R);
                end
            else
                if(right) %#ok
                    thetaT = pi/2;
                else
                    thetaT = pi/2;
                end
            end

            if(~insideBottom)
                if(right)
                    thetaB = asin((y2Min - OD2)/R);
                else
                    thetaB = pi - asin((y2Min - OD2)/R);
                end
            else
                if(right)
                    thetaB = -pi/2;
                else
                    thetaB = 3*pi/2;
                end
            end
    
            tCircle = linspace(thetaB, thetaT, 100);
            
            plot(R*cos(tCircle) + OD1, R*sin(tCircle) + OD2, ...
                '-r', 'HandleVisibility','off');

            % fix the axes
            axis equal
            box on
            set(gca,'XTick',[])
            set(gca,'YTick',[])            
            xlim(xLim);
            ylim(yLim);
            legend('Location','southeast');
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
            geom.Origin = [1;1];
            geom.OriginDisk = [-4,0];

            % geom.y2Min = -3;  % full in the slit
            % geom.y2Max = 2;
            
            % geom.y2Min = -0.5; % wider than the slit
            % geom.y2Max = -0.2;

            % geom.y2Min = -3; % cuts top
            % geom.y2Max = -0.5;

            geom.y2Min = 0.5; % cuts bottom
            geom.y2Max = 3;

            geom.L1 = 3;
            % geom.LeftRight = 'Left';
            geom.LeftRight = 'Right';

        end

    end



end