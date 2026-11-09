classdef HalfSpaceQuad < SpectralSpectral
    % The upper half space ([-inf,inf] x [y2Min,inf] in Cartesian coordinates.
    % Note that this is the same as HalfSpace except we use a quadratic map
    % in the first coordinate, rather than a quotient map.
    % Use HalfSpace unless you have a good reason not to.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L1 - a mapping parameter such that half of the y1 points lie in
    %  [-L1 + y10, L1 + y10]
    %  L2 - a mapping parameter such that half of the y2 points lie in
    %  [y2Min, y2Min + L2]
    %  (y10 = 0) - the 'origin' in y1

    properties
        y2Min = 0; % position of wall in y2
        L1  % mapping parameter in y1
        L2  % mapping parameter in y2
        y10 = 0; % 'origin' in y2
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = HalfSpaceQuad(Geometry)

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1    = Geometry.L1;
            this.L2    = Geometry.L2;

            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end
            if(isfield(Geometry,'y2Min'))
                this.y2Min = Geometry.y2Min;
            end

            % set that it's Cartesian coordinates
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
            % quadratic map from [-1,1] to [-inf,inf]
            % shifted by y10
            [y1,dy1,dx,ddx,dddx,ddddx] = QuadMap(x1,this.L1,inf);
            y1 = y1 + this.y10;
        end

        function xf = CompSpace1(this,y1)
            % inverse quadratic map
            xf  = InvQuadMap(y1-this.y10,this.L1,inf);
        end

        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % quotient map from [-1,1] to [y2Min,inf]
            [th,dth,dx,ddx,dddx,ddddx] = QuotientMap(x2,this.L2,...
                this.y2Min,inf);
        end

        function x2 = CompSpace2(this,y2)
            % inverse quotient map
            x2  = InvQuotientMap(y2,this.L2,this.y2Min,inf);
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
            O1 = this.y10;
            O2 = this.y2Min;

            % set sensible axis limits
            xLim = [-1.1 * 3*L1 + O1, 1.1 * 3*L1 + O1];
            yLim = [O2 - 0.1*3*L2, 1.1 * 3*L2 + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;

            % plot indicator lines for L1, L2
            hold on
            plot([-L1 + O1, -L1 + O1],yLim,'--m','DisplayName', ' \pm Li + Origin(i)');
            plot([L1 + O1,L1 + O1],yLim,'--m','HandleVisibility','off');
            plot(xLim,[L2 + O2, L2 + O2],'--m','HandleVisibility','off');
            
            % plot indicator for wall
            plot(xLim,[O2,O2],'-r','DisplayName', 'y2Min');

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

            geom.N = [12;20];
            geom.L1 = 3;
            geom.L2 = 5;
            geom.y10 = -2;
            geom.y2Min = 3;

        end

    end

end