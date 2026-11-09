classdef BoxTanY2 < SpectralSpectral
    % A box with corners defined by (y1Min,y2Min), (y1Max,y2Min),
    % (y1Min,Y2Max), and (y1Max,y2Max).
    % Alternatively, specify the lengths of the two sides, L1 and L2, and
    % (optionally) the Origin, which is the bottom left of the box
    %
    %
    % This differs from BOX by using a tangent-based map in y2 to cluster
    % the points around a given value - essentially an adapted grid.
    % If alpha1>>1, points cluster around alpha2.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  (y1Min = 0)      - left edge
    %  y1Max (or L1)    - right edge
    %  (y2Min = 0)      - bottom
    %  y2Max (or L2)    - top
    %  (Origin = [0,0]) - bottom left of the box when using L1, L2
    %  alpha1  - parameter for tan map
    %  alpha2  - parameter for tan map


    properties
        y1Min = 0
        y2Min = 0
        y1Max
        y2Max

        alpha1 % map parameter
        alpha2 % map parameter
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = BoxTanY2(Geometry)

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

            % map parameters
            this.alpha1 = Geometry.alpha1;
            this.alpha2 = Geometry.alpha2;

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
            %
            % This is the same as in Box

            disp('Computing Convolution matrices...');

            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            if(isfield(shapeParams,'shape'))
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

                % set up the subshape and intersect it with the box

                geom = ShiftGeometryOrigin(shapeParams,[y1(i);y2(i)]);

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
        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % tan map from [-1,1] to [y2Min,y2Max]
            [y,dy,dx,ddx,dddx,ddddx] = TanIntervalMap(x2,this.alpha1,this.alpha2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)
            % inverse tan map
            % if alpha1>>1, points cluster around alpha2.
            x2  = InvTanIntervalMap(y2,this.alpha1,this.alpha2,this.y2Min,this.y2Max);
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
            alpha2 = this.alpha2;

            % set sensible axis limits
            xLim = [y1Min + O1 - 0.1*L1, y1Max + O1 + 0.1*L1];
            yLim = [y2Min + O2 - 0.1*L2, y2Max + O2 + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            plot(xLim,[alpha2,alpha2],'--m','DisplayName','alpha2');

            % plot the box
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
            geom.alpha1 = 20;
            geom.alpha2 = 1.5;

        end

    end
    
end