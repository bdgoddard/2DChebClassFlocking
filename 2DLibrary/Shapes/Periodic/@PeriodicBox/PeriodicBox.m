classdef PeriodicBox < SpectralFourier
    % A (half) periodic box with corners defined by
    % (y1Min,y2Min), (y1Max,y2Min), (y1Min,Y2Max), and (y1Max,y2Max).
    % Alternatively, specify the lengths of the two sides, L1 and L2, and
    % (optionally) the Origin, which is the bottom left of the box
    % The box is periodic in the second coordinate.
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
        y2Min = 0 % bottom (first periodic boundary)
        y1Max     % right
        y2Max     % top (second periodic boundary)
        L2        % length of periodic dimension
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = PeriodicBox(Geometry)
            this@SpectralFourier(Geometry.N(1),Geometry.N(2));

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

            % define the periodic length
            this.L2 = this.y2Max - this.y2Min;

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

            int = ComputeIntegrationVector@SpectralFourier(this);

            area = ( this.y1Max - this.y1Min ) * ( this.y2Max - this.y2Min );

            % display the error in the area if requested
            if(nargin == 2 && displayArea)

                if(area == 0)
                    disp(['PeriodicBox: Area is zero, Absolute error is: ',...
                        num2str(area-sum(this.Int))]);
                else
                    disp(['PeriodicBox: Error of integration of area (ratio): ',...
                        num2str(1-sum(this.Int)/area)]);
                end

            end

        end

        function M_conv = ComputeConvolutionMatrix_Pointwise(this,f,shapeParams)
            % Convolution \int_D(x) f(x-y)g(y) dy
            % Here D(x) is the intersection of the subshape defined in
            % shapeParams with the box

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

                % can modify disc by specifying the origin
                % need to shift y1/2/Min/Max for box
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
                    fP = f(GetDistance(this,dy1,dy2));
                else
                    fP = f(GetDistance1(this,dy1),GetDistance2(this,dy2));
                end

                for iF = 1:fDim
                    M_conv(i,:,iF) = (Int_i.*fP(:,iF)')*IP_i;
                end

            end

            M_conv(isnan(M_conv)) = 0;

            disp('...done.');

            close(wb);

        end % convolution

        
        function M_IntConv = ComputeIntegrationConvolutionMatrix(this,fInt,fConv,shapeParams)
            % Computes the matrix that does a weighted integral in the
            % first coordinate (with function fInt) and a convolution in
            % the second coordinate (with fConv).  The convolution can be
            % performed on a subdomain specified by shapeParams.

            if(nargin<4)
                shapeParams = struct;
            end

            % Weighted integration in y1 (spectral)

            % set up spectral line
            geom1.N = this.N1;
            geom1.yMin = this.y1Min;
            geom1.yMax = this.y1Max;
            sLine = SpectralLine(geom1);

            % compute integration weights
            Int1 = sLine.ComputeIntegrationVector;
            ws = sLine.Pts.y;
            gLine = fInt(ws);
            I1Line = Int1.* gLine';
            I1Rep = repmat(I1Line,[this.N1,1]);
            I1 = sparse(kron(I1Rep,eye(this.N2)));

            % Convolution in y2 (Fourier)

            % set up Fourier line
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);

            % compute convolution (which requires the integration vector)
            fLine.ComputeIntegrationVector;
            % format: C2 = fLine.ComputeConvolutionMatrix(this,f,shapeParams,parent)
            C2Line = fLine.ComputeConvolutionMatrix(fConv,shapeParams);
            C2 = sparse(kron(eye(this.N1),C2Line));

            % Full integral
            M_IntConv = C2*I1;
        end


        function [M_ConvConv,C1,C2,C1Line,C2Line] = ComputeConvolutionConvolutionMatrix(this,f1,f2,shapeParams1,shapeParams2)
            % Computes the matrix that does convolutions in both
            % coordinates (with f1 and f2, respectively).  Both can be
            % performed on subdomains specified by shapeParams1 and
            % shapeParams2.

            % Convolution in y1 (spectral) with cutoff

            % set up spectral line
            geom1.N = this.N1;
            geom1.yMin = this.y1Min;
            geom1.yMax = this.y1Max;
            sLine = SpectralLine(geom1);

            % compute convolution in y1
            sLine.ComputeIntegrationVector;
            C1Line = sLine.ComputeConvolutionMatrix_Pointwise(f1,shapeParams1);
            C1 = sparse(kron(C1Line,eye(this.N2)));

            % Convolution in y2 (Fourier)

            % set up Fourier line
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);

            % compute convolution (which requires the integration vector)
            fLine.ComputeIntegrationVector;
            C2Line = fLine.ComputeConvolutionMatrix(f2,shapeParams2);
            C2 = sparse(kron(eye(this.N1),C2Line));

            % Full convolution
            M_ConvConv = C1*C2;
        end


        function [Int1,Int2] = Compute1DIntegrationMatrices(this)
            % Partial integration in each of the two coordinates

            % Spectral integration
            geom1.N = this.N1;
            geom1.yMin = this.y1Min;
            geom1.yMax = this.y1Max;
            sLine = SpectralLine(geom1);
            I1 = sLine.ComputeIntegrationVector;
            Int1 = sparse(kron(I1,eye(this.N2)));

            % Fourier integration
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            fLine = FourierLine(geom2);
            I2 = fLine.ComputeIntegrationVector;
            Int2 = sparse(kron(eye(this.N1),I2));
        end

        function [M_conv_1D,M_conv] = ComputeConvolutionMatrix_1D_1(this,f,shapeParams)
            % Convolution in the first coordinate with kernel f and on a
            % subdomain defined by shapeParams.

            if(nargin<2)
                shapeParams = struct;
            end

            % Create line that corresponds to the y1 line
            geom1.N = this.N1;
            geom1.yMin = this.y1Min;
            geom1.yMax = this.y1Max;
            y1Line = SpectralLine(geom1);

            % Compute 1D convolution
            y1Line.ComputeIntegrationVector;

            % this will default to y1Line.ComputeConvolutionMatrix if
            % there's no shapeParams
            M_conv_1D = y1Line.ComputeConvolutionMatrix_Pointwise(f,shapeParams);

            % Copy for each y2
            M_conv = sparse(kron(M_conv_1D,eye(this.N2)));

        end

        function [M_conv_1D,M_conv] = ComputeConvolutionMatrix_1D_2(this,f,shapeParams)
            % Convolution in the second coordinate with kernel f and on a
            % subdomain defined by shapeParams.

            if(nargin<2)
                shapeParams = struct;
            end

            % Create line that corresponds to the y2 line
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            y2Line = FourierLine(geom2);

            % Compute 1D convolution
            y2Line.ComputeIntegrationVector;
            % no pointwise option as it's not needed for periodic
            M_conv_1D = y2Line.ComputeConvolutionMatrix(f,shapeParams);

            % Copy for each y1
            M_conv = sparse(kron(eye(this.N1),M_conv_1D));

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
            % linear map from [0,1] to [y2Min,y2Max]
            % note that this is Fourier so defined on [0,1]
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(x2,this.y2Min,this.y2Max);
        end

        function x2 = CompSpace2(this,y2)
            % inverse linear map
            x2 = InvLinearMap01(y2,this.y2Min,this.y2Max);
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
            plot([y1Min,y1Min],[y2Min,y2Max],'-r','DisplayName', 'Non-periodic boundary');
            plot([y1Max,y1Max],[y2Min,y2Max],'-r','HandleVisibility','off');
            plot([y1Min,y1Max],[y2Min,y2Min],'--r','DisplayName', 'Periodic boundary');
            plot([y1Min,y1Max],[y2Max,y2Max],'--r','HandleVisibility','off');

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