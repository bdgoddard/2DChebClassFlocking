classdef Wedge < Polar_SpectralSpectral
    % A wedge of a disc of radius R, between th1 and th2, centred at
    % Origin.
    % Can also be used to create a wedge of an annulus by specifying R_in,
    % the inner radius.
    %
    % Inputs:
    %   N - [N1,N2] number of points in the two coordinates
    %   R (or R_out) - (outer) radius of the disc
    %   th1, th2 - the start and end angles
    %   (Origin = [0,0]) - origin of the circle
    %
    % Alternatively, specify:
    %   y1Min - R_in
    %   y1Max - R_out
    %   y2Min - th1
    %   y2Max - th2

    properties
        R_in = 0;  % inner radius
        R_out      % outer radius
        th1,th2    % angle range
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = Wedge(Geometry)

            this@Polar_SpectralSpectral(Geometry.N(1),Geometry.N(2));


            if(isfield(Geometry,'y1Min'))
                % set through specifying all four options in the usual
                % Shape format
                this.R_in  = Geometry.y1Min;
                this.R_out = Geometry.y1Max;
                this.th1   = Geometry.y2Min;
                this.th2   = Geometry.y2Max;
            else
                % specify R and theta
                if(isfield(Geometry,'R'))
                    this.R_out    = Geometry.R;
                else
                    this.R_out    = Geometry.R_out;
                end

                if(isfield(Geometry,'R_in'))
                    this.R_in     = Geometry.R_in;
                end

                this.th1      = Geometry.th1;
                % only need the mod-2pi difference between the angles
                this.th2      = Geometry.th1 + mod(Geometry.th2 - Geometry.th1,2*pi);
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % set the origin
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

        end
    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,t1Odd,t2Odd,displayArea)

            % compute the integration weights
            if(nargin == 1)
                int  = ComputeIntegrationVector@Polar_SpectralSpectral(this);
            elseif(nargin == 2)
                int  = ComputeIntegrationVector@Polar_SpectralSpectral(this,t1Odd);
            elseif(nargin == 3)
                int  = ComputeIntegrationVector@Polar_SpectralSpectral(this,t1Odd,t2Odd);
            end

            % the exact area
            area = (this.th2-this.th1)/2*(this.R_out^2 - this.R_in^2);

            if(nargin == 4 && displayArea)
                if(area == 0)
                    disp(['Wedge: Area is zero. Error of integration of area(abs): ', num2str(sum(this.Int))]);
                else
                    disp(['Wedge: Error of integration of area(ratio): ', num2str(1-sum(this.Int)/area)]);
                end
            end

        end

        function Ind  = ComputeIndices(this)
            % compute the indices for boundaries

            Ind      = GetIndicesWedge(this);
            this.Ind = Ind;

        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------
        
        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % Chebyshev linear map for radius
            [r,dr,dx,ddx,dddx,ddddx] = LinearMap(x1,this.R_in,this.R_out);
        end

        function x1 = CompSpace1(this,r)
            % Inverse linear map
            x1  = InvLinearMap(r,this.R_in,this.R_out);
        end

        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % Chebyshev linear map for angle (note that this is not a
            % periodic construction)
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.th1,this.th2);
        end

        function x2 = CompSpace2(this,y2)
            % Inverse linear map
            x2 = InvLinearMap(y2,this.th1,this.th2);
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
            R_in = this.R_in; %#ok<*PROPLC>
            R_out = this.R_out;
            th1 = this.th1;
            th2 = this.th2;
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * R_out + O1, 1.1 * R_out + O1];
            yLim = [-1.1 * R_out + O2, 1.1 * R_out + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot indicator circles
            tCircle = linspace(0, 2*pi, 100);
            plot(R_out * cos(tCircle) + O1, R_out * sin(tCircle) + O2, ...
                 '--k', 'DisplayName', 'Circles of radius R_{in}, R_{out}');
            plot(R_in * cos(tCircle) + O1, R_in * sin(tCircle) + O2, ...
                 '--k', 'HandleVisibility','off');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicator arrows
            % arrow for R_out
            p1 = [O1, O2];
            theta = (th1+th2)/2;
            p2 = [O1 + R_out * cos(theta), O2 + R_out * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2 - 0.02 * R_out, 'R out', ...
                'FontSize',20,'Color','b')
            % arrow for R_in
            if(R_in > 0)
                theta = (th1+th2)/3;
                p2 = [O1 + R_in * cos(theta), O2 + R_in * sin(theta)];
                dp = p2 - p1;
                plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
                text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2 + 0.02 * R_out, 'R in', ...
                    'FontSize',20,'Color','b')
            end

            % plot indicators for the two angles
            tth1 = linspace(0, mod(th1,2*pi), 100);
            tth2 = linspace(0, mod(th2,2*pi), 100);
            R1 = R_in + 0.2*(R_out-R_in);
            R2 = R_in + 0.3*(R_out-R_in);
            plot(R1 * cos(tth1) + O1, R1 * sin(tth1) + O2, ...
                '-m', 'DisplayName', 'th1, th2');
            plot(R2 * cos(tth2) + O1, R2 * sin(tth2) + O2, ...
                '-m', 'HandleVisibility','off')
            % line at theta = 0
            plot([O1,xLim(2)],[O2,O2],':m','HandleVisibility','off')


            % plot the boundary
            tCircle = linspace(th1, th2, 100);
            plot(R_out*cos(tCircle) + O1, R_out*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');
            if(R_in > 0)
                plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                    '-r', 'HandleVisibility','off');
            end
            plot([R_in*cos(th1) + O1, R_out*cos(th1) + O1],[R_in*sin(th1) + O2, R_out*sin(th1) + O2], ...
                    '-r', 'HandleVisibility','off');
            plot([R_in*cos(th2) + O1, R_out*cos(th2) + O1],[R_in*sin(th2) + O2, R_out*sin(th2) + O2], ...
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

            geom.N = [10;12];
            geom.R_in = 0.2;
            geom.R_out = 1;
            geom.th1 = pi/6;
            geom.th2 = 1.5*pi;
            geom.Origin = [-0.2,0.5];

        end

    end

end