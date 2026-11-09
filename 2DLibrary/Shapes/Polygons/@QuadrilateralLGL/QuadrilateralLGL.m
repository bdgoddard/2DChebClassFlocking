classdef QuadrilateralLGL < M1SpectralSpectralLGL
    % A _convex_ quadrilateral specified by four corners, y1, y2, y3, y4.
    % Note: this is exactly the same as QUADRILATERAL except we use
    % Legende-Gauss-Lobatto nodes and weights.  If in doubt, use
    % QUADRILATERAL.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  Y = [y1;y2;y3;y4] - with yi being 1x2 vectors of the Cartesian
    %  coordinates of the corners


    properties
        Y % corner points, mapped in order from (-1,-1), (1,-1), (1,1), (-1,1)
        % Y(i,:) = (y_i1,y_i2)
        A % mapping matrix from (0,0), (1,0), (1,1), (0,1) -> Y

    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = QuadrilateralLGL(Geometry)

            % ensure that the points are in an appropriate cyclic order
            Geometry.Y = SortQuadrilateralIndices(Geometry.Y);

            % note that this uses LGL rather than Chebyshev (it's worse)
            this@M1SpectralSpectralLGL(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.Y   = Geometry.Y;
            this.A   = GetQuadrilateralParameters(this.Y);

            % we're in Cartesian coordinates
            this.polar = 'cart';
            
            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
        end

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function [int,area] = ComputeIntegrationVector(this,displayArea)
            
            int  = ComputeIntegrationVector@M1SpectralSpectral(this);
            
            % exact area
            x = this.Y(:,1);
            y = this.Y(:,2);
            area = 1/2*abs( ( x(1)*y(2) - y(1)*x(2) ) + ( x(2)*y(3) - y(2)*x(3) ) ...
                + ( x(3)*y(4) - y(3)*x(4) ) + ( x(4)*y(1) - y(4)*x(1) ));

            % display the error in the area if requested
            if(nargin == 2 && displayArea)
                if(area == 0)
                    disp(['Quadrilateral: Error of integration of area(=0): ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['Quadrilateral: Error of integration of area(ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end
            end
        end

        function Ind    = ComputeIndices(this)
            % find the indices
            Ind      = GetIndicesQuadrilateral(this);
            this.Ind = Ind;
        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            % see e.g. https://www.particleincell.com/2012/quad-interpolation/

            % shorthands
            n  = length(x1);
            O  = ones(n,1);

            % shift [-1,1] to [0,1] via linear map
            x1 = (x1+1)/2;
            x2 = (x2+1)/2;

            % Jacobians
            Dx1 = 1/2;
            Dx2 = 1/2;

            % apply mapping from [0,1] x [0,1] to quadrilateral determined
            % by corners Y
            h = this.A*[O';x1';x2';(x1.*x2)'];
            y1_kv = h(1,:)';  y2_kv = h(2,:)';

            % construct Jacobian and Hessians
            if(nargout >= 3)
                % Jacobian
                J        = zeros(n,2,2);
                J(:,1,1) = (this.A(1,2)*O + this.A(1,4)*x2) * Dx1;
                J(:,1,2) = (this.A(1,3)*O + this.A(1,4)*x1) * Dx2;
                J(:,2,1) = (this.A(2,2)*O + this.A(2,4)*x2) * Dx1;
                J(:,2,2) = (this.A(2,3)*O + this.A(2,4)*x1) * Dx2;
            end

            if(nargout >= 4)
                % Hessian of y1 map
                dH1        = zeros(n,2,2);
                dH1(:,1,2) = this.A(1,4)*O * Dx1 * Dx2;
                dH1(:,2,1) = this.A(1,4)*O * Dx1 * Dx2;
            end

            if(nargout >= 4)
                % Hessian of y2 map
                dH2        = zeros(n,2,2);
                dH2(:,1,2) = this.A(2,4)*O * Dx1 * Dx2;
                dH2(:,2,1) = this.A(2,4)*O * Dx1 * Dx2;
            end

        end

        function [x1,x2] = CompSpace(this,y1,y2)

            % map onto [0,1] x [0,1]
            alpha = this.A(1,:);
            beta = this.A(2,:);

            % a x2^2 + b x2 + c = 0
            a = alpha(4)*beta(3) - alpha(3)*beta(4);
            b = alpha(4)*beta(1) - alpha(1)*beta(4) + alpha(2)*beta(3) ...
                - alpha(3)*beta(2) + y1*beta(4) - y2*alpha(4);
            c = alpha(2)*beta(1) - alpha(1)*beta(2) + y1*beta(2) - y2*alpha(2);

            if(abs(a)<10^-10) % should be a=0, but include a tolerance
                w2 = -c./b;
            else
                det = sqrt(b.^2-4*a.*c);
                w2 = (-b + det)/(2*a);
            end
            w1 = (y1 - alpha(1) - alpha(3)*w2)./(alpha(2) + alpha(4)*w2);

            % map onto [-1,1] x [-1,1]
            x1 = 2*w1 - 1;
            x2 = 2*w2 - 1;
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
            Y = this.Y; %#ok<*PROPLC>

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

            % plot the boundary
            plot([Y(1,1),Y(2,1)],[Y(1,2),Y(2,2)],'-r','DisplayName','Boundary');
            plot([Y(2,1),Y(3,1)],[Y(2,2),Y(3,2)],'-r','HandleVisibility','off');
            plot([Y(3,1),Y(4,1)],[Y(3,2),Y(4,2)],'-r','HandleVisibility','off');
            plot([Y(4,1),Y(1,1)],[Y(4,2),Y(1,2)],'-r','HandleVisibility','off');

            % plot the corners
            plot(Y(1,1),Y(1,2),'or','MarkerSize',5,'MarkerFaceColor','r','DisplayName','y1');
            plot(Y(2,1),Y(2,2),'ob','MarkerSize',5,'MarkerFaceColor','b','DisplayName','y2');
            plot(Y(3,1),Y(3,2),'om','MarkerSize',5,'MarkerFaceColor','m','DisplayName','y3');
            plot(Y(4,1),Y(4,2),'oc','MarkerSize',5,'MarkerFaceColor','c','DisplayName','y4');

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
            geom.Y = [0, -1;
                      1, 1;
                      -2, 0.5;
                      4, -2];

        end

    end
    
end

