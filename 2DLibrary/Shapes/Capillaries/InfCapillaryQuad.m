classdef InfCapillaryQuad < InfCapillaryGeneral
    % An infinite (in the first coordinate) capillary, which is finite in
    % the second coordinate between y2Min and y2Max.
    % Note this only differs from InfCapillary through the use of a
    % quadratic map (rather than square root) in the first coordinate.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L1 - mapping parameter for the first coordinate
    %  y2Min - bottom
    %  y2Max - top
    %  y10 - origin in first coordinate
    %  (Conv = false) - true/false to set up for FMT-like convolution

    % copied from InfCapillaryGeneral for ease of reference
    % properties
    %     L1      % mapping parameter for infinite dimension
    %     y2Min   % bottom
    %     y2Max   % top
    %     L1Conv  % convolution mapping parameter
    %     L2Conv  % convolution mapping parameter
    %     y10 = 0 % origin in first cooordinate
    % end


    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = InfCapillaryQuad(Geometry)
            this@InfCapillaryGeneral(Geometry);
        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % quadratic map from [-1,1] to [-inf,inf] centred at y10
            [y1,dy1,dx,ddx,dddx,ddddx] = QuadMap(x1,this.L1,inf);
            y1 =  this.y10 + y1;
        end

        function xf = CompSpace1(this,y1)
            % inverse quadratic map
            xf  = InvQuadMap(y1 - this.y10,this.L1,inf);
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
            y2Min = this.y2Min;
            y2Max = this.y2Max;
            L2 = y2Max - y2Min;
            O1 = this.y10;

            % set sensible axis limits
            xLim = [-1.1 * 3*L1 + O1, 1.1 * 3*L1 + O1];
            yLim = [y2Min - 0.1*L2, y2Max + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot domain limits
            plot(xLim,[y2Min,y2Min],'-r','DisplayName', 'Walls');
            plot(xLim,[y2Max,y2Max],'-r','HandleVisibility','off');

            % plot indicator lines for L1
            plot([-L1 + O1, -L1 + O1],[y2Min,y2Max],'--m','DisplayName', 'Origin(1) \pm L1');
            plot([L1 + O1,L1 + O1],[y2Min,y2Max],'--m','HandleVisibility','off');

            % plot the origin
            plot(O1,0,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','Origin')

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

            geom.N = [12;20];
            geom.y2Min = -1;
            geom.y2Max = 3;
            geom.L1 = 3;
            geom.y10 = -2;

        end

    end


end