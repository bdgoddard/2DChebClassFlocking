classdef Disc < Polar_SpectralFourier
    % A disc of radius R, centred at Origin.
    % Note that this is defined in polar coordinates with y1 the radius (on
    % (0,R]) and y2 the angle (on [0,2pi)).
    % 
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  R - the radius of the disc
    %  (Origin = [0,0]) - [O1,O2] the origin of the disc
    %  (spherical = false) - whether to do spherically-weighted integration

    properties
        R  % radius
        sphere = false;  % whether to do 'spherical' integration
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        
        function this = Disc(Geometry)

            this@Polar_SpectralFourier(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.R      = Geometry.R;

            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;
            end

            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)

            if(this.sphere)
                % include the spherical weighting
                int = ComputeIntegrationVector@Polar_SpectralFourier(this,true);
                int = 2*int.*sqrt(this.R^2 - this.Pts.y1_kv.^2)';

                % check a common error
                if(sum(imag(int)) > 0)
                    error('Disc:ComputeIntegrationVector:Imaginary part in integration vector');
                end

                % assign to Int
                this.Int = int;

                % the exact area
                area     = 4/3*pi*this.R^3;
            else

                int = ComputeIntegrationVector@Polar_SpectralFourier(this);
                
                % the exact area
                area = pi*this.R^2;
            end

            % display the error in the area if requested
            if(nargin == 2 && displayArea)
                disp(['Disc: Error of integration of area (ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);
            end

        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            % standard linear map from [-1,1] to [-R,R]
            % note that Polar_SpectralFourier will cut this down to only
            % the positive radius values, after having doubled N1
            % so the points are singly defined, but there is more effort
            % needed to compute the differentiation and integration
            [r,dr,dx,ddx,dddx,ddddx] = LinearMap(x,-this.R,this.R);
        end

        function xf = CompSpace1(this,r)
            % inverse linear map
            xf = InvLinearMap(r,-this.R,this.R);
        end

        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(~,xT)
            % linear map from [0,1] to [0,2pi]
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end

        function xf = CompSpace2(~,th)
            % inverse linear map
            xf = InvLinearMap01(th,0,2*pi);
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

            % set sensible axis limits
            xLim = [-1.4 * R + O1, 1.4 * R + O1];
            yLim = [-1.4 * R + O2, 1.4 * R + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','Origin')

            % plot indicator arrows
            % arrow for R
            p1 = [O1, O2];
            theta = pi/6;
            p2 = [O1 + R * cos(theta), O2 + R * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.03 * R, p1(2) + dp(2)/2 - 0.03 * R, 'R', ...
                'FontSize',20,'Color','b')

            % plot the boundary
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'DisplayName', 'Boundary');


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
            geom.R = 1;
            geom.Origin = [0.2,0.7];

        end

    end


end