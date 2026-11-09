classdef InfCapillarySkewed < InfCapillary
    % An infinite (in the first coordinate) capillary, which is finite in
    % the second coordinate between y2Min and y2Max.  The coordinate grid
    % is skewed at an angle alpha.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  L1 - mapping parameter for the first coordinate
    %  y2Min - bottom
    %  y2Max - top
    %  y10 - origin in first coordinate
    %  (Conv = false) - true/false to set up for FMT-like convolution
    %  (alpha = pi/2) - skew of grid

    % copied from InfCapillaryGeneral for ease of reference
    % properties
    %     L1      % mapping parameter for infinite dimension
    %     y2Min   % bottom
    %     y2Max   % top
    %     L1Conv  % convolution mapping parameter
    %     L2Conv  % convolution mapping parameter
    %     y10 = 0 % origin in first cooordinate
    % end

    properties
        alpha = pi/2;  % skew of grid
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = InfCapillarySkewed(Geometry)

            % set parameters after skewing
            Geometry.y2Min = Geometry.y2Min/sin(Geometry.alpha);
            Geometry.y2Max = Geometry.y2Max/sin(Geometry.alpha);

            this@InfCapillary(Geometry);

            % set parameter
            this.alpha  = Geometry.alpha;

            % ensure that we are in a skewed Cartesian grid
            this.polar  = 'cartSkewed';

            % set up points using method from Shape (and PhysSpace below)
            InitializationPts(this);
        end

    end

    methods (Access = public)

        %------------------------------------------------------------------
        % Overridden initialisations
        %------------------------------------------------------------------

        function ptsCart = GetCartPts(this,pts_y1,pts_y2)

            % allow using a structure for Pts
            if(nargin == 1)
                pts_y1 = this.Pts.y1_kv;
                pts_y2 = this.Pts.y2_kv;
            end

            % get the skewed grid
            [ptsCart.y1_kv,ptsCart.y2_kv] = SkewedGrid(pts_y1,pts_y2,this.alpha);
        end


        function pts = GetInvCartPts(this,ptsCart_y1,ptsCart_y2)
        
            % should just return the points as we're not in polar
            if(strcmp(this.polar,'cartSkewed'))
                [pts.y1_kv,pts.y2_kv] = InvSkewedGrid(ptsCart_y1,ptsCart_y2,this.alpha);
            else
                exc = MException('HalfSpaceSkewed:GetInvCartPts','select cartSkewed');
                throw(exc);
            end

        end

        function [Int,Int1,Int2]    = ComputeIntegrationVector(this)
            %This gives the integration in the real (cartesian) space!
            [Int,Int1,Int2] = ComputeIntegrationVector@SpectralSpectral(this);
            Int             = Int*sin(this.alpha);
            Int2            = Int2*sin(this.alpha);
        end

        function Diff = ComputeDifferentiationMatrix(this)
            % include the skew in the differentiation

            DiffSk = ComputeDifferentiationMatrix@SpectralSpectral(this);

            ATInv = zeros(2,2);
            ATInv(1,1) = 1;
            ATInv(2,1) = -1/tan(this.alpha);
            ATInv(2,2) = 1/sin(this.alpha);
            %A(1,2) = cos(this.alpha);
            %A(2,2) = sin(this.alpha);

            % y = A*yt
            %nabla_y f = (A^(T))^(-1)* nabla_y'*f

            h         = zeros(this.M,this.M,2);
            h(:,:,1)  = DiffSk.Dy1;
            h(:,:,2)  = DiffSk.Dy2;

            hh        = MultScalarMatrVecOp(ATInv,h);
            Diff.Dy1  = hh(:,:,1);
            Diff.Dy2  = hh(:,:,2);

            %Diff.Dy1  = ATInv(1,1)*DiffSk.Dy1 + ATInv(1,2)*DiffSk.Dy2;
            %Diff.Dy2  = ATInv(2,1)*DiffSk.Dy1 + ATInv(2,2)*DiffSk.Dy2;

            Diff.grad = sparse([Diff.Dy1;Diff.Dy2]);

            H          = zeros(this.M,this.M,2,2);
            H(:,:,1,1) = DiffSk.DDy1;
            H(:,:,1,2) = DiffSk.Dy1Dy2;
            H(:,:,2,1) = DiffSk.Dy1Dy2;
            H(:,:,2,2) = DiffSk.DDy2;

            w          = MultScalarMatrOp(ATInv,H);

            ht         = w(:,:,1,2);
            w(:,:,1,2) = w(:,:,2,1);
            w(:,:,2,1) = ht;

            %     w          = permute(w,[1 2 4 3]);
            w          = MultScalarMatrOp(ATInv,w);

            ht         = w(:,:,1,2);
            w(:,:,1,2) = w(:,:,2,1);
            w(:,:,2,1) = ht;

            Diff.DDy1    = w(:,:,1,1);
            Diff.Dy1Dy2  = w(:,:,1,2);
            Diff.DDy2    = w(:,:,2,2);%DiffSk.DDy2;%
            %
            %Diff.DDy1 = ;

        end

        function M_conv = ComputeConvolutionMatrix(this,f,saveBool,ptsCheck)
            % allow pointwise convolution
            if(nargin>=4)
                M_conv  = ComputeConvolutionMatrix@InfCapillary(this,f,ptsCheck);
            else
                M_conv  = ComputeConvolutionMatrix@InfCapillary(this,f);
            end

            M_conv = M_conv*sin(this.alpha);

            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end

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
            alpha = this.alpha;
            y2Min = this.GetCartPts(0,this.y2Min).y2_kv;
            y2Max = this.GetCartPts(0,this.y2Max).y2_kv;
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
            plot([-L1 + O1, (-L1 + O1) + L2/tan(alpha)],[y2Min,y2Max],'--m','DisplayName', 'L1');
            plot([L1 + O1, (L1 + O1) + L2/tan(alpha)],[y2Min,y2Max],'--m','HandleVisibility','off');

            % plot indicator for alpha
            talpha = linspace(0, mod(alpha,2*pi), 100);
            R = L2/6;
            plot(R * cos(talpha) + L1 + O1, R * sin(talpha) + y2Min, ...
                '-b', 'DisplayName', 'alpha');

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

            % get the defaults from InfCapillary
            geom = InfCapillary.DocumentationPlotDefaultGeometry();
            
            % set a default angle
            geom.alpha = pi/3;

        end

    end

end
