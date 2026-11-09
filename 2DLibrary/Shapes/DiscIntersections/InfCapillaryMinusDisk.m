classdef InfCapillaryMinusDisk < handle
    % An infinite (in the first coordinate) capillary with walls at y1Min
    % and y2Max, with a disc centred at Origin removed. Note that the
    % Origin must lie inside the capillary, and the disc must not
    % completely disect the capillary.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radius of the disc to be removed
    %  L1 - mapping parameter for the first coordinate
    %  y2Min, y2Max - positions of the upper and lower walls
    %  (Origin = [0,0]) - origin of the disc

    properties
        SubShape % stores the subshapes
        Int % integration vector
        y2Min, y2Max  % wall positions
        R % disc radius
        Origin % disc origin
        L1 % mapping parameter
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = InfCapillaryMinusDisk(Geometry)

            % set parameters
            this.L1 = Geometry.L1;
            this.y2Min = Geometry.y2Min;
            this.y2Max = Geometry.y2Max;
            this.R = Geometry.R;
            this.Origin = Geometry.Origin;

            % construct on a case-by-case basis

            if(Geometry.y2Max - Geometry.y2Min <= 2*Geometry.R)
                % capillary is narrower than the disc
                error('Not yet implemented');

            elseif(((Geometry.Origin(2)-Geometry.y2Min) >= Geometry.R) && ...
                    ((Geometry.y2Max - Geometry.Origin(2)) >= Geometry.R))
                % disc fits into the capillary
                Geometry.TopBottom = 'Top';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{1} = StripMinusDisk(Geometry);
                Geometry.y2Wall    = Geometry.y2Min;
                Geometry.TopBottom = 'Bottom';
                this.SubShape{2} = StripMinusDisk(Geometry);

            elseif(Geometry.Origin(2) == Geometry.y2Min)
                % disc is centred on the bottom wall
                Geometry.TopBottom = 'Top';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{1} = StripMinusDisk(Geometry);

            elseif(((Geometry.Origin(2)-Geometry.y2Min) < Geometry.R) ...
                    && ((Geometry.Origin(2)-Geometry.y2Min) >= 0))
                % disc is within R of the bottom wall inside the capillary
                Geometry.TopBottom = 'Top';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{1} = StripMinusDisk(Geometry);

                Geometry.TopBottom = 'Bottom';
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                Geometry.y2Wall    = Geometry.y2Min;
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);

                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);

            elseif(Geometry.Origin(2) == Geometry.y2Max)
                % disc is centred on the top wall
                Geometry.TopBottom = 'Bottom';
                Geometry.y2Wall    = Geometry.y2Min;
                this.SubShape{1} = StripMinusDisk(Geometry);

            elseif((abs(Geometry.Origin(2)-Geometry.y2Max) < Geometry.R) ...
                    && ((Geometry.Origin(2)-Geometry.y2Max) <= 0))
                % disc is within R of the top wall inside the capillary
                Geometry.TopBottom = 'Bottom';
                Geometry.y2Wall    = Geometry.y2Min;
                this.SubShape{1} = StripMinusDisk(Geometry);

                Geometry.TopBottom = 'Top';
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);

                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);
            else
                error('Not yet implemented');
            end

        end

        %------------------------------------------------------------------
        % Shape-like functionality
        %------------------------------------------------------------------

        function ptsCart = GetPts(this)
            % concatenate the points
            ptsCart.y1_kv = [];
            ptsCart.y2_kv = [];
            for i = 1:length(this.SubShape)
                ptsCart.y1_kv = [ptsCart.y1_kv;this.SubShape{i}.Pts.y1_kv];
                ptsCart.y2_kv = [ptsCart.y2_kv;this.SubShape{i}.Pts.y2_kv];
            end
        end

        function ptsCart = GetCartPts(this)
            % use the Subshapes to get the Cartesian points and concatenate
            % them
            ptsCart.y1_kv = [];
            ptsCart.y2_kv = [];

            for i = 1:length(this.SubShape)
                ptsCartSub    = this.SubShape{i}.GetCartPts();
                ptsCart.y1_kv = [ptsCart.y1_kv;ptsCartSub.y1_kv];
                ptsCart.y2_kv = [ptsCart.y2_kv;ptsCartSub.y2_kv];
            end
        end

        function int = ComputeIntegrationVector(this)
            % concatenate the integration weights
            int = [];
            for i = 1:length(this.SubShape)
                int = [int,this.SubShape{i}.ComputeIntegrationVector()]; %#ok<AGROW>
            end
        end

        function PlotGridLines(this)
            %  plot gridlines for each SubShape
            for i = 1:length(this.SubShape)
                this.SubShape{i}.PlotGridLines(); hold on;
            end
        end

        function PlotGrid(this)
            % plot grid points for each SubShape
            for i = 1:length(this.SubShape)
                this.SubShape{i}.PlotGrid(); hold on;
            end
        end

        function plot(this,V,opts)
            % plot using the subshapes
            if(nargin==2)
                opts = [];
            end
            mOld = 1;
            for i = 1:length(this.SubShape)
                this.SubShape{i}.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
                m = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2 - 1;
                this.SubShape{i}.plot(V(mOld:m),opts);
                mOld = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2;
                hold on;
            end
            hold on;
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
            O1 = this.Origin(1);
            O2 = this.Origin(2);
            
            centre = ( (O2-y2Min) >= R) && ((y2Max - O2) >= R);
            top = ( abs(y2Max - O2) < R );
            bottom = ( abs(y2Min - O2) < R );

            % set sensible axis limits
            xLim = [-1.1 * 3*L1 +  O1, 1.1 * 3*L1 + O1];
            yLim = [y2Min - 0.5*L2, y2Max + 0.5*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator line for L1
            plot([O1+L1,O1+L1],[y2Min,y2Max],'--m','DisplayName','Origin(1) \pm L1');
            plot([O1-L1,O1-L1],[y2Min,y2Max],'--m','HandleVisibility','off');

            % plot circle
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName','Disc and walls');
            % plot walls
            plot(xLim,[y2Max,y2Max],'--k', 'HandleVisibility','off');
            plot(xLim,[y2Min,y2Min],'--k', 'HandleVisibility','off');

            % plot domain decomposition
            if(top || bottom)
                plot(xLim,[O2,O2],':k','DisplayName','Domain split');
            end

            % plot boundary
            % horizontal lines
            if(centre)
                plot(xLim,[y2Min,y2Min],'-r', 'DisplayName', 'Boundary');
                plot(xLim,[y2Max,y2Max],'-r', 'HandleVisibility','off');
                tCircle = linspace(0, 2*pi, 100);
            end
            if(top)
                theta = 3*pi/2 + acos((y2Max - O2)/R);
                plot(xLim,[y2Min,y2Min],'-r', 'DisplayName', 'Boundary');
                plot([-R*sin(theta),xLim(2)],[y2Max,y2Max],'-r', 'HandleVisibility','off');
                plot([xLim(1), R*sin(theta) + O1],[y2Max,y2Max],'-r', 'HandleVisibility','off');
                tCircle = linspace(theta+pi, theta + 5*pi/2, 100);
            end
            if(bottom)
                theta = 3*pi/2 + acos((O2-y2Min)/R);
                plot(xLim,[y2Max,y2Max],'-r', 'DisplayName', 'Boundary');
                plot([-R*sin(theta) + O1,xLim(2)],[y2Min,y2Min],'-r', 'HandleVisibility','off');
                plot([xLim(1),  R*sin(theta) + O1],[y2Min,y2Min],'-r', 'HandleVisibility','off');
                tCircle = linspace(theta + 3*pi/2, theta, 100);
            end
            % arc
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '-r', 'HandleVisibility','off');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','OriginDisk')

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
            
            % geom.y2Min = -0.7; % bottom
            % geom.y2Max = 3;

            geom.y2Min = -3; % top
            geom.y2Max = 0.7;

            % geom.y2Min = -3; % centre
            % geom.y2Max = 3;

            geom.L1 = 3;
            geom.Origin = [0,0];

        end

    end


end

