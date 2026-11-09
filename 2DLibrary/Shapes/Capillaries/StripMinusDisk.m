classdef StripMinusDisk < Polar_M1SpectralSpectral
    % An infinite, horizontal strip with a disc of radius R centred at
    % Origin removed.  The other edge of the strip is at y2Wall.  The half
    % disc must fit in the strip, otherwise use HalfStripMinusDisk.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radius of the disc to be removed
    %  L1 - mapping parameter for radial coordinate; half the points lie in
    %       [R,R+L1]
    %  y2Wall - position of lower wall
    %  (TopBottom = 'Bottom') - 'Top' or 'Bottom' for whether strip is at
    %        the Top or Bottom of the disc.
    %  (Origin = [0,0]) - origin of the disc

    properties
        R  % radius of disc
        y2Wall % distance of the wall from the origin
        L1  % mapping parameter in infinite coordinate
        TopBottom = 'Bottom' % location of wall
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = StripMinusDisk(Geometry)

            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.y2Wall     = Geometry.y2Wall;

            if(isfield(Geometry,'Origin'))
                this.Origin     = Geometry.Origin;
            end

            if(isfield(Geometry,'TopBottom'))
                this.TopBottom = Geometry.TopBottom;
            end

            % check that constructing the shape makes sense
            err = false;
            top = strcmp(this.TopBottom,'Top');
            if(top && (this.y2Wall < this.Origin(2)))
                err = true;
            end
            if(~top && (this.y2Wall > this.Origin(2)))
                err = true;
            end
            if(err)
                exc = MException('StripMinusDisk:Constructor', ...
                    'y2Wall and TopBottom must be compatible.');
                throw(exc);
            end

            if(abs(this.Origin(2)-this.y2Wall) <= this.R)
                exc = MException('StripMinusDisk:Constructor', ...
                    'Half disc must fit in the Strip.  Use HalfStripMinusDisk.');
                throw(exc);
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

            O = ones(size(x1));

            % polar coordinate angle depending on orientation
            if(strcmp(this.TopBottom,'Top'))
                [y2_kv,dy2] = LinearMap(x2,0,pi);
            elseif(strcmp(this.TopBottom,'Bottom'))
                [y2_kv,dy2] = LinearMap(x2,pi,2*pi);
            end

            % radial distance from the origin for each y2
            rd  = abs((this.Origin(2)-this.y2Wall)./sin(y2_kv));
            rd((y2_kv==pi)|(y2_kv==2*pi)) =  inf;

            % do some robustness scaling
            if((this.Origin(2)-this.y2Wall-this.R) > 3*this.L1)
                L1_r = this.L1*ones(size(rd));
            else
                L1_r          = this.L1*(rd-this.R)./(this.L1*3+rd-this.R);
                L1_r(rd==inf) = this.L1;
            end

            % quotient map for the radiuse
            [y1_kv,dy1]   = QuotientMap(x1,L1_r,O*this.R,rd);

            % Jacobian and Hessians
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
            exc = MException('StripMinusDisk:CompSpace','not yet implemented');
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
            O1 = this.Origin(1);
            O2 = this.Origin(2);
            L2 = abs(O2 - y2Wall);
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
            plot(xLim,[y2Wall,y2Wall],'-r', 'DisplayName', 'Border');
            plot([xLim(1),O1-R],[O2,O2],'-r','HandleVisibility','off');
            plot([O1+R,xLim(2)],[O2,O2],'-r','HandleVisibility','off');

            % plot curved edge
            if(top)
                tCircle = linspace(0, pi, 100);
            else
                tCircle = linspace(pi, 2*pi, 100);
            end

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
            geom.L1 = 1.5;
            geom.Origin = [-1,2];

            % geom.y2Wall = 0;  % option 1
            % geom.TopBottom = 'Bottom';

            geom.y2Wall = 4;  % option 2
            geom.TopBottom = 'Top';

        end

    end


end