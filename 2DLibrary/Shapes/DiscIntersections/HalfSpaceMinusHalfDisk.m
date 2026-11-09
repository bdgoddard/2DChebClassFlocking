classdef HalfSpaceMinusHalfDisk < Polar_M1SpectralSpectral
    % An half space from Origin(2) in the second coordinate with half a
    % disc of radius R and centred at Origin removed.  Which half space is
    % determined by TopBottom.  This is in polar coordinates.
    % Can also specify Rmax, which results in half an annulus.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radius of the disc to be removed
    %  (Rmax = inf) - outer radius
    %  L1 - mapping parameter for radial coordinate; half the points lie in
    %       [R,R+L1]
    %  (Origin = [0,0]) - origin of the disc

    properties
        R   % radius of the disc
        L1  % radial mapping parameter
        Rmax = Inf  % outer radius
        TopBottom = 'Top'  % 'Top' or 'Bottom' half space
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = HalfSpaceMinusHalfDisk(Geometry)

            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;

            if(isfield(Geometry,'Origin'))
                this.Origin     = Geometry.Origin;
            end
            if(isfield(Geometry,'Rmax'))
                this.Rmax       = Geometry.Rmax;
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

            % quotient map from R to Rmax
            [y1_kv,dy1] = QuotientMap(x1,this.L1,this.R,this.Rmax);

            % linear map for angle
            if(strcmp(this.TopBottom,'Top'))
                [y2_kv,dy2] = LinearMap(x2,0,pi);
            elseif(strcmp(this.TopBottom,'Bottom'))
                [y2_kv,dy2] = LinearMap(x2,pi,2*pi);
            end

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
            exc = MException('HalfSpaceMinusDisk:CompSpace','not yet implemented');
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
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * 3*L1 + O1, 1.1 * 3*L1 + O1];
            yLim = [O2 - 0.1*3*L1, 1.1 * 3*L1 + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator line for L1+R
            tCircle = linspace(0, pi, 100);
            plot((L1+R)*cos(tCircle) + O1, (L1+R)*sin(tCircle) + O2, ...
                    '--m', 'DisplayName','L1+R');
            
            % plot circle
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName','Wall and Disc');
            p1 = [O1, O2];
            theta = 3*pi/4;
            p2 = [O1 + R * cos(theta), O2 + R * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * 3*L1, p1(2) + dp(2)/2 + 0.02 * 3*L1, 'R', ...
                'FontSize',20,'Color','b')

            % plot wall
            plot(xLim,[O2,O2],'-k', 'HandleVisibility','off');


            % plot wall
            plot([R + O1,xLim(2)],[O2,O2],'-r', 'DisplayName', 'Boundary');
            plot([xLim(1),O1 - R],[O2,O2],'-r', 'HandleVisibility','off');
            % plot half disc
            tCircle = linspace(0, pi, 100);
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
            geom.L1 = 3;
            geom.TopBottom = 'Top';
            geom.Origin = [-1,2];

        end

    end

end