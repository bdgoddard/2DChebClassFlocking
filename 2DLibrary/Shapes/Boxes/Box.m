classdef Box < SpectralSpectral
    % A box with corners defined by (y1Min,y2Min), (y1Max,y2Min),
    % (y1Min,Y2Max), and (y1Max,y2Max).
    % Alternatively, specify the lengths of the two sides, L1 and L2, and
    % (optionally) the Origin, which is the bottom left of the box
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

        function this = Box(Geometry)

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

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

        function M_conv = ComputeConvolutionMatrix_Pointwise(this,f,shapeParams)
            % Convolution \int_D(x) f(x-y)g(y) dy
            % Here D(x) is the intersection of the subshape defined in
            % shapeParams with the box

            disp('Computing Convolution matrices...');

            % determine if the kernel uses a distance or the two
            % coordinates
            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            if(isfield(shapeParams,'shape'))
                % determine the shape we're using
                subshapeFn = str2func(shapeParams.shape);

                % Disc, which will need to be decomposed into segments,
                % triangles and boxes; NT is the number of points for the
                % triangle
                if(strcmp(shapeParams.shape,'Disc'))
                    if(~isfield(shapeParams,'NT'))
                        shapeParams.NT = shapeParams.N;
                    end
                end

            else
                % whole box
                M_conv = ComputeConvolutionMatrix(this,f);
                return;

            end

            Pts = this.Pts;

            if(useDistance)
                fPTemp = f(GetDistance(this,Pts.y1_kv,Pts.y2_kv));
            else
                fPTemp = f(Pts.y1_kv,Pts.y2_kv);
            end

            fDim = size(fPTemp,2);

            N1  = this.N1;  N2  = this.N2;
            N1N2 = N1*N2;
            M_conv = zeros(N1*N2,N1*N2,fDim);

            wb = waitbar(0,'Computing Convolution Matrices');

            y1 = Pts.y1_kv;
            y2 = Pts.y2_kv;

            for i=1:N1N2

                waitbar(i/N1N2,wb);

                % shift the origin to the current point
                geom = ShiftGeometryOrigin(shapeParams,[y1(i);y2(i)]);

                % set up the subshape and intersect it with the box
                subshapeFull = subshapeFn(geom);
                subShape = Intersect(this,subshapeFull);

                % determine the points in the intersected domain
                subConvShapePts = subShape.pts;

                % the interpolation onto D(x) from the full domain
                IP_i     = SubShapePts(this,subConvShapePts);

                % the integration weights for the subshape
                Int_i    = subShape.int;

                % f(x-y)
                dy1 = y1(i) - subShape.pts.y1_kv;
                dy2 = y2(i) - subShape.pts.y2_kv;

                if(useDistance)
                    fP          = f(GetDistance(this,dy1,dy2));
                else
                    fP          = f(dy1,dy2);
                end

                for iF = 1:fDim
                    M_conv(i,:,iF) = (Int_i.*fP(:,iF)')*IP_i;
                end

            end

            M_conv(isnan(M_conv)) = 0;

            disp('...done.');

            close(wb);

        end % convolution

        function [M_conv,M_conv_1D] = ComputeConvolutionMatrix_Pointwise_1D_1(this,f,shapeParams)
            % convolution only in the first dimension

            % Create line that corresponds to the y1 line
            geom1.N = this.N1;
            geom1.yMin = this.y1Min;
            geom1.yMax = this.y1Max;
            y1Line = SpectralLine(geom1);

            % Compute 1D convolution
            y1Line.ComputeIntegrationVector;
            M_conv_1D = y1Line.ComputeConvolutionMatrix_Pointwise(f,shapeParams);

            % Copy for each y2
            M_conv = sparse(kron(M_conv_1D,eye(this.N2)));

        end


        function [int,area] = ComputeIntegrationVector(this,displayArea)

            int = ComputeIntegrationVector@SpectralSpectral(this);

            area = ( this.y1Max - this.y1Min ) * ( this.y2Max - this.y2Min );

            % display the error in the area if requested
            if(nargin == 2 && displayArea)

                if(area == 0)
                    disp(['Box: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['Box: Error of integration of area (ratio): ',...
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

        function [y2,dy2,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % linear map from [-1,1] to [y2Min,y2Max]
            [y2,dy2,dx,ddx,dddx,ddddx] = LinearMap(x2,this.y2Min,this.y2Max);
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

            % set sensible axis limits
            xLim = [y1Min - 0.1*L1, y1Max + 0.1*L1];
            yLim = [y2Min - 0.1*L2, y2Max + 0.1*L2];

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