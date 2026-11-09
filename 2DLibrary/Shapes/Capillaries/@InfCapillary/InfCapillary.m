classdef InfCapillary < InfCapillaryGeneral & ConvolutionFiniteSupport
    % An infinite (in the first coordinate) capillary, which is finite in
    % the second coordinate between y2Min and y2Max.
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

        function this = InfCapillary(Geometry)
            this@InfCapillaryGeneral(Geometry);

            % set up the FMT convolution
            %There are four numberings:
            % (here, there is only one domain,simplifying the numbering, as
            % all identifiers = 1)
            %(1) the numberings of the domains. The points with identifyer
            %     X are given by mark_id(:,X)
            this.mark_id = true(this.M,1);
            %(2) the numbering of the domains in the first dimension. Points
            %in the first dimension belonging to the part with identifyer Y
            %are given by mark_id_1(:,Y)
            this.mark_id_1 = true(this.N1,1);
            %(3) the numbering of the domains in the second dimension. Points
            %in the second dimension belonging to the part with identifyer Z
            %are given by mark_id_2(:,Z)
            this.mark_id_2 = true(this.N2,1);
            %mark12 makes a link between identifyers of dimensions 1,2 and
            % the global identifyer. Points which are in the first dimension
            % located in a domain with identifyer Y and in the second
            % direction in a domain with identifyer Z, belong to the global
            % domain with the identifyer X = mark_12(Y,Z).
            this.mark_12   = ones(1,1);
        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        % function for FMT convolutions
        AD = ComputeConvolutionFiniteSupport(this,area,weights,pts,params)

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------
 
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % square root map from [-1,1] to [-inf,inf] centred at y10
            [y1,dy1,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);
            y1 =  this.y10 + y1;
        end

        function xf = CompSpace1(this,y1)
            % inverse square root map
            xf  = InvSqrtMap(y1 - this.y10,this.L1,inf);
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
            plot([-L1 + O1, -L1 + O1],[y2Min,y2Max],'--m','DisplayName', 'L1');
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