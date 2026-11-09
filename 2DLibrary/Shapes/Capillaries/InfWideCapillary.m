% File: InfWideCapillary.m
% Created Date: 16 Oct 2025
% Author: Jake Skelton
% Date Modified: Thu Oct 16 2025
% Modification of InfCapillaryGeneral for a capillary wide in the finite direction

classdef InfWideCapillary < SpectralSpectral
    % An infinite (in the first coordinate) capillary, which is finite in
    % the second coordinate between y2Min and y2Max.
    % This differs from InfCapillary in that it uses QuotientMap in the
    % finite direction rather than LinearMap.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L1 - mapping parameter for the first coordinate
    %  L2 - mapping parameter for the first coordinate
    %  y2Min - bottom
    %  y2Max - top
    %  y10 - origin in first coordinate

    properties
        L1 (1,1) double % Length scale in infinite direction
        L2 (1,1) double % Length scale in finite direction
        y2Min (1,1) double % bottom
        y2Max (1,1) double % top
        y10 (1,1) double % y1 origin
    end
    
    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = InfWideCapillary(Geometry)
          
            arguments
                Geometry.N (1,2) double {mustBePositive,mustBeInteger}
                Geometry.L1 (1,1) double {mustBePositive}
                Geometry.L2 (1,1) double {mustBePositive}
                Geometry.y2Min (1,1) double
                Geometry.y2Max (1,1) double
                Geometry.y10 (1,1) double = 0
            end

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            % set parameters
            this.L1      = Geometry.L1;
            this.L2      = Geometry.L2;
            this.y2Min   = Geometry.y2Min;
            this.y2Max   = Geometry.y2Max;
            this.y10     = Geometry.y10;
            
            % set that we're in Cartesian coordinates
            this.polar = 'cart';
            
            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
            
        end
    end
    
    methods (Access = public)

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
        
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % quotient map from [-1,1] to [y2Min,y2Max]
            [th,dth,dx,ddx,dddx,ddddx] = QuotientMap(...
                x2, this.L2, this.y2Min, this.y2Max);
        end
        
        function x2 = CompSpace2(this,y2)
            % inverse quotient map
            x2 = InvQuotientMap(y2,this.L2,this.y2Min,this.y2Max);
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
            L2 = this.L2;
            y2Min = this.y2Min;
            y2Max = this.y2Max;
            Ly2 = y2Max - y2Min;
            O1 = this.y10;

            % set sensible axis limits
            xLim = [-1.1 * 4*L1 + O1, 1.1 * 4*L1 + O1];
            yLim = [y2Min - 0.1*Ly2, y2Max + 0.1*Ly2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator lines for L1
            plot([-L1 + O1, -L1 + O1],[y2Min,y2Max],'--m','DisplayName', 'Origin(1) \pm L1');
            plot([L1 + O1,L1 + O1],[y2Min,y2Max],'--m','HandleVisibility','off');

            % plot indicator line for L2
            plot(xLim,[L2 + y2Min,L2 + y2Min],'--m','DisplayName','y2Min + L2');

            % plot the origin
            plot(O1,0,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','Origin')

            % plot boundary
            plot(xLim,[y2Min,y2Min],'-r','DisplayName', 'Boundary');
            plot(xLim,[y2Max,y2Max],'-r','HandleVisibility','off');

            % fix the axes
            axis equal
            box on
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            xlim(xLim);
            ylim(yLim);
            legend('Location','west');
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
            geom.L1 = 4;
            geom.L2 = 0.5;
            geom.y10 = -2;

            geom = namedargs2cell(geom);

        end

    end


end