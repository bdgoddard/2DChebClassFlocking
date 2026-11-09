classdef BoxFD < FiniteDifferenceFiniteDifference
    % A box with corners defined by (y1Min,y2Min), (y1Max,y2Min),
    % (y1Min,Y2Max), and (y1Max,y2Max).
    % Alternatively, specify the lengths of the two sides, L1 and L2, and
    % (optionally) the Origin, which is the bottom left of the box.
    % NOTE: THIS USES FINITE DIFFERENCE
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  (y1Min = 0)      - left edge
    %  y1Max (or L1)    - right edge
    %  (y2Min = 0)      - bottom
    %  y2Max (or L2)    - top
    %  (Origin = [0,0]) - bottom left of the box when using L1, L2


    properties
        y1Min = 0 % left
        y2Min = 0 % bottom
        y1Max     % right
        y2Max     % top
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = BoxFD(Geometry)

            this@FiniteDifferenceFiniteDifference(Geometry.N(1),Geometry.N(2));

            % set parameters
            if(isfield(Geometry,'L1'))
                % defined by the lengths of the two sides
                
                % set left, right, top, and bottom
                if(isfield(Geometry,'Origin'))
                    this.y1Min = Geometry.Origin(1);
                    this.y2Min = Geometry.Origin(2);
                end
                this.y1Max      = this.y1Min + Geometry.L1;
                this.y2Max      = this.y2Min + Geometry.L2;

            else
                % left, right, top, and bottom are given

                this.y1Min = Geometry.y1Min;
                this.y1Max = Geometry.y1Max;
                this.y2Min = Geometry.y2Min;
                this.y2Max = Geometry.y2Max;
            end

            % set that we're in Cartesian coordinates
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

            int = ComputeIntegrationVector@FiniteDifferenceFiniteDifference(this);

            area = ( this.y1Max - this.y1Min ) * ( this.y2Max - this.y2Min );

            % display the error in the area if requested
            if(nargin == 2 && displayArea)

                if(area == 0)
                    disp(['BoxFD: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['BoxFD: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end

            end

        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % linear map from [-1,1] to [y1Min,y1Max]
            [y1,dy1,dx,ddx,dddx,ddddx] = LinearMap(x1,this.y1Min,this.y1Max);
        end

        function xf = CompSpace1(this,y1)
            % inverse linear map
            xf  = InvLinearMap(y1,this.y1Min,this.y1Max);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % linear map from [-1,1] to [y2Min,y2Max]
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)
            % inverse linear map
            x2  = InvLinearMap(y2,this.y2Min,this.y2Max);
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
            y1Min = this.y1Min; %#ok<*PROPLC>
            y1Max = this.y1Max;
            y2Min = this.y2Min;
            y2Max = this.y2Max;
            L1 = y1Max - y1Min;
            L2 = y2Max - y2Min;
            O1 = this.Origin(1);
            O2 = this.Origin(2);

            % set sensible axis limits
            xLim = [y1Min + O1 - 0.1*L1, y1Max + O1 + 0.1*L1];
            yLim = [y2Min + O2 - 0.1*L2, y2Max + O2 + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;

            % plot the box
            hold on
            plot([y1Min,y1Min],[y2Min,y2Max],'-r','DisplayName', 'Boundary');
            plot([y1Max,y1Max],[y2Min,y2Max],'-r','HandleVisibility','off');
            plot([y1Min,y1Max],[y2Min,y2Min],'-r','HandleVisibility','off');
            plot([y1Min,y1Max],[y2Max,y2Max],'-r','HandleVisibility','off');

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
            geom.y1Min = -1;
            geom.y1Max = 3;
            geom.y2Min = 0;
            geom.y2Max = 2;

        end

    end

end