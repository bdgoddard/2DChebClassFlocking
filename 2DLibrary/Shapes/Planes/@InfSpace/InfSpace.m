classdef InfSpace < SpectralSpectral
    % The whole space in Cartesian coordinates.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L1 - a mapping parameter such that half of the y1 points lie in
    %  [-L1 + y10, L1 + y10]
    %  L2 - a mapping parameter such that half of the y2 points lie in
    %  [-L2 + y20, L2 + y20]
    %  (y10 = 0) - the 'origin' in y1
    %  (y20 = 0) - the 'origin' in y2

    properties
        AD = [];  % used for FMT calculations
        L1 % mapping parameter
        L2 % mapping parameter
        y10 = 0; % 'origin' in y1
        y20 = 0; % 'origin' in y2
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------
        
        function this = InfSpace(Geometry)
            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1    = Geometry.L1;
            this.L2    = Geometry.L2;

            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end
            if(isfield(Geometry,'y20'))
                this.y20 = Geometry.y20;
            end
            
            % set that it's Cartesian coordinates
            this.polar = 'cart';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
        end
    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams,parent)
            %{
            Strategy:
            We wish to compute \int_D(x) f(x-y) g(y) dy, for some domain D
            (e.g. disc, annulus) centred at x.
            Let z = y-x then dy = dz, y = x + z and 
            y \in D(x) => z + x \in D(x) or z \in D(0) and we get
            \int_D(0) f(-z) g(x+z) dz
            
            In matlab this is
            M_conv = Int_D(0) .* f(-z) .* IP_D(x)
            %}

            if((nargin>=4) && parent)
                M_conv = ComputeConvolutionMatrix@M1SpectralSpectral(this,f,[]);
                return;
            end

            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            disp('Computing Convolution matrices...');

            n1 = this.N1;   n2 = this.N2;

            if(isfield(shapeParams,'R'))
                disp('Disc');
                subShape = Disc(shapeParams);
                rRange = (-1:0.05:1)';
            elseif(isfield(shapeParams,'RMin'))
                if(isfield(shapeParams,'RMax'))
                    disp('Annulus');
                    subShape = Annulus(shapeParams);
                    rRange = (-1:0.05:1)';
                else
                    disp('InfAnnulus');
                    subShape = InfAnnulus(shapeParams);
                    rRange = (-1:0.05:0.9)';
                end
            else
                disp('InfDisc');
                subShape = InfDisc(shapeParams);
                rRange = (-0.9:0.05:0.9)';
            end

            % Integration and interpolation on D(0)
            Int_i    = subShape.ComputeIntegrationVector();
            subShape.ComputeInterpolationMatrix(rRange,(0:0.02:1)',true,true);

            % f(-z)
            subShapeCartPts = Pol2CartPts(subShape.Pts);
            invSubShapePts = invertPts(subShapeCartPts,'cart');

            if(useDistance)
                fP       = f(this.GetDistance(invSubShapePts.y1_kv,invSubShapePts.y2_kv));
            else
                fP       = f(invSubShapePts.y1_kv,invSubShapePts.y2_kv);
            end

            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end

            nf = size(fP,2);

            IntRep = Int_i(ones(1,nf),:);

            M_conv = zeros(n1*n2,n1*n2,nf);

            figure
            subShape.plot(fP);
            title('Interpolation of Convolution function'); drawnow;


            y1_kv=this.Pts.y1_kv;
            y2_kv=this.Pts.y2_kv;

            hw = waitbar(0,'Computing Convolution matrices...');

            for i=1:(n1*n2)

                waitbar(i/(n1*n2),hw);

                Y10 = y1_kv(i);
                Y20 = y2_kv(i);

                % Compute points in D(x) = (Y10,Y20))
                subConvShapePts = shiftPointsCart(subShapeCartPts,Y10,Y20);
                % and the interpolation onto D(x) from the full domain
                IP          = SubShapePts(this,subConvShapePts);

                IntFTemp = (IntRep.*fP.')*IP;
                M_conv(i,:,:) = IntFTemp.';

            end

            close(hw);

            M_conv(isnan(M_conv)) = 0;

            if(nDim>2)
                nfTemp = size(fP_old);
                nfTemp = nfTemp(2:end);
                M_conv = reshape(M_conv,[ n1*n2, n1*n2, nfTemp ]);
            end

        end


        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            % square root map from [-1,1] to [-inf,inf]
            % shifted by y10
            [y1,dy1,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);
            y1 = y1 + this.y10;
        end

        function x1 = CompSpace1(this,y1)
            % inverse square root map
            x1  = InvSqrtMap(y1-this.y10,this.L1,inf);
        end

        function [y2,dy2,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)
            % square root map from [-1,1] to [-inf,inf]
            % shifted by y20
            [y2,dy2,dx,ddx,dddx,ddddx] = SqrtMap(x2,this.L2,inf);
            y2 = y2 + this.y20;
        end

        function x2 = CompSpace2(this,y2)
            % inverse square root map
            x2  = InvSqrtMap(y2-this.y20,this.L2,inf);
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
            O2 = this.y20;

            % set sensible axis limits
            xLim = [-1.1 * 4*L1 + O1, 1.1 * 4*L1 + O1];
            yLim = [-1.1 * 3*L2 + O2, 1.1 * 3*L2 + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;

            % plot indicator lines for L1, L2
            hold on
            plot([-L1 + O1, -L1 + O1],yLim,'--m','DisplayName', 'Origin(i) \pm Li');
            plot([L1 + O1,L1 + O1],yLim,'--m','HandleVisibility','off');
            plot(xLim,[-L2 + O2, -L2 + O2],'--m','HandleVisibility','off');
            plot(xLim,[L2 + O2, L2 + O2],'--m','HandleVisibility','off');

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',10,'MarkerFaceColor','k','DisplayName','Origin')

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

            geom.N = [12;20];
            geom.L1 = 4;
            geom.L2 = 3;
            geom.y10 = -2;
            geom.y20 = 3;

        end

    end


    methods (Access = private)
        % test functions

        ComputeConvolutionMatrix_Test(this);
        testAll(this);

    end
end