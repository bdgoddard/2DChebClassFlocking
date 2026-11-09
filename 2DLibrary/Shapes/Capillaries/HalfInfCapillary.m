classdef HalfInfCapillary < SpectralSpectral
    % A capillary that is half infinite in the first coordinate and finite
    % in the second, with bounds given by y2Min and y2Max.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L1 - mapping parameter for first coordinate.  Half the points lie in
    %       [0,L1]
    %  y2Min - bottom
    %  y2Max - top

    properties
        L1    % mapping parameter
        y2Min % bottom
        y2Max % top
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = HalfInfCapillary(Geometry)

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1      = Geometry.L1;
            this.y2Min   = Geometry.y2Min;
            this.y2Max   = Geometry.y2Max;

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

        function [int,area] = ComputeIntegrationVector(this)

            % standard integration vector
            int  = ComputeIntegrationVector@SpectralSpectral(this);
            
            % provide the exact area
            area = Inf;
        end


        function M_conv = ComputeConvolutionMatrix(this,f,saveBool,ptsCheck)

            % allow to specify pointwise convolution
            if(nargin>=4)
                M_conv  = ComputeConvolutionMatrix_Pointwise(this,f,ptsCheck);
            else
                M_conv  = ComputeConvolutionMatrix_Pointwise(this,f);
            end

            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
        end
        
        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------
        
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % quotient map from [-1,1] to [0,inf]
            [y1,dy1,dx,ddx,dddx,ddddx] = QuotientMap(x1,this.L1,0);
        end

        function xf = CompSpace1(this,y1)
            % inverse quotient map
            xf  = InvQuotientMap(y1,this.L1,0);
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
            y2Min = this.y2Min; %#ok<*PROPLC>
            y2Max = this.y2Max;
            L1 = this.L1;
            L2 = y2Max - y2Min;

            % set sensible axis limits
            xLim = [- 0.1*L1, 3.1*L1];
            yLim = [y2Min - 0.1*L2, y2Max + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator line for L1
            plot([L1, L1],[y2Min,y2Max],'--m','DisplayName', 'L1');

            % plot the boundary
            plot([0,0],[y2Min,y2Max],'-r','DisplayName', 'Boundary');
            plot([0,xLim(2)],[y2Min,y2Min],'-r','HandleVisibility','off');
            plot([0,xLim(2)],[y2Max,y2Max],'-r','HandleVisibility','off');

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
            geom.L1 = 3;
            geom.y2Min = 0;
            geom.y2Max = 2;

        end

    end


end