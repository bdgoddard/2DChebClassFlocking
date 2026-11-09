classdef Triangle < M1SpectralSpectral
    % A triangle specified by the three corners, y1, y2, y3.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  Y = [y1;y2;y3] - with yi being 1x2 vectors of the Cartesian
    %  coordinates of the corners

    properties
        Y % the corner points
        A % matrix used in PhysSpace
    end
    
    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = Triangle(Geometry)

            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.Y   = Geometry.Y;
            this.A   = GetTriangleParameters(this.Y);

            % we're in Cartesian coordinates
            this.polar = 'cart';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);

            % map is not invertible since one side is mapped to a corner
            this.invertibleMap = false; 

        end

    end

    methods (Access = public)
                
        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)

            % use the standard integration vector
            int  = ComputeIntegrationVector@M1SpectralSpectral(this);

            % explicit area
            y = this.Y;
            area = 1/2*abs((y(1,1)-y(3,1)) * (y(2,2)-y(1,2))-...
                           (y(1,1)-y(2,1)) * (y(3,2)-y(1,2)));       

            % display the error in the area if requested
            if(nargin == 2 && displayArea)
                if(area == 0)
                    disp(['Triangle: Error of integration of area(=0): ',...
                                        num2str(area-sum(this.Int))]);
                else
                    disp(['Triangle: Error of integration of area(ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);
                end
            end
        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)  

            % this is a standard map from a box to a triangle
            n  = length(x1);
            O  = ones(n,1);
            A = this.A; %#ok<*PROPLC>

            h = A*[O';x1';x2';(x1.*x2)'];
            y1_kv = h(1,:)';  y2_kv = h(2,:)';

            % compute the Jacobian and Hessians
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = A(1,2)*O + A(1,4)*x2;
                J(:,1,2) = A(1,3)*O + A(1,4)*x1;
                J(:,2,1) = A(2,2)*O + A(2,4)*x2;
                J(:,2,2) = A(2,3)*O + A(2,4)*x1;
            end

            if(nargout >= 4)
                dH1        = zeros(n,2,2);     
                dH1(:,1,2) = A(1,4)*O;
                dH1(:,2,1) = A(1,4)*O;
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);            
                dH2(:,1,2) = A(1,4)*O;
                dH2(:,2,1) = A(1,4)*O;
            end

        end

        function [x1,x2] = CompSpace(~,~)  %#ok<STOUT>
             % we don't implent the computational space map
             exc = MException('Triangle:CompSpace','not yet implemented');
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

            % shorthand variables
            Y = this.Y;

            % set sensible axis limits
            xMin = min(Y(:,1));
            xMax = max(Y(:,1));
            yMin = min(Y(:,2));
            yMax = max(Y(:,2));
            % add in some minimal padding
            xPad = 0.1 * max(abs(xMin),1);
            yPad = 0.1 * max(abs(xMin),1);
            xLim = [xMin - xPad, xMax + yPad];
            yLim = [yMin - xPad, yMax + yPad];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot the corners
            hold on
            plot(Y(1,1),Y(1,2),'oc','MarkerSize',5,'MarkerFaceColor','c','DisplayName','y1');
            plot(Y(2,1),Y(2,2),'ob','MarkerSize',5,'MarkerFaceColor','b','DisplayName','y2');
            plot(Y(3,1),Y(3,2),'om','MarkerSize',5,'MarkerFaceColor','m','DisplayName','y3');

            % plot the boundary
            plot([Y(1,1),Y(2,1)],[Y(1,2),Y(2,2)],'-r','DisplayName','Boundary');
            plot([Y(2,1),Y(3,1)],[Y(2,2),Y(3,2)],'-r','HandleVisibility','off');
            plot([Y(3,1),Y(1,1)],[Y(3,2),Y(1,2)],'-r','HandleVisibility','off');

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
            geom.Y = [0, 0;
                1, 1;
                -2, 0.5];

        end

    end


end

