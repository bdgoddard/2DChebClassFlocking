classdef FluctuatingChannel < M1SpectralSpectral
    % A finite channel, of length L1, where the top and bottom are 
    % defined by functions, ft and fb.
    %
    % Inputs:
    %  N - [N1,N2] the number of points in the two coordinates
    %  ft - function name for top; input is (x,params), output is y
    %  fb - function name for bottom; input is (x,params), output is y
    %  ftparams - parameters for ft
    %  fbparams - parameters for fb
    
    properties
        L1  % mapping parameter for infinite domain

        ft  % function to give top
        fb  % function to give bottom
        ftParams  % parameters for ft
        fbParams  % parameters for fb

    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = FluctuatingChannel(Geometry)

            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            % set parameters
            this.L1 = Geometry.L1;

            % set top and bottom functions and parameters
            this.ft = str2func(Geometry.ft);
            this.fb = str2func(Geometry.fb);
            this.ftParams = Geometry.ftParams;
            this.fbParams = Geometry.fbParams;

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

        function Ind    = ComputeIndices(this)
            % custom indices
            Ind      = GetIndicesFluctuatingChannel(this);
            this.Ind = Ind;
        end

        %------------------------------------------------------------------
        % Maps
        %------------------------------------------------------------------

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            
            % shorthand parameters
            L1      = this.L1; %#ok<*PROPLC>
            ft      = this.ft;
            fb      = this.fb;
            ftParams = this.ftParams;
            fbParams = this.fbParams;
            n  = length(x1);

            % map is given by
            % y1 = (x1+1)*L1/2;
            % y2 = bottom(x1) + (x2+1).*(top(x1)-bottom(x1))/2;

            % map first coordinate from [-1,1] to [0,L1]
            [y1_kv,dy1dx1] =  LinearMap(x1,0,L1);

            % find top and bottom as functions of x
            [top, dftdx1, dftdxx1] = ft(x1,ftParams);
            [bottom, dfbdx1, dfbdxx1] = fb(x1,fbParams);

            % figure
            % x1
            % plot(x1,top,'r');
            % hold on
            % plot(x1,bottom,'b');

            % pointwise linear map from bottom to top
            [y2_kv,dy2dx2] =  LinearMap(x2,bottom,top);

            % determine Jacobian and Hessians
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = dy1dx1;
                %J(:,1,2) = zeros(n,1);
                J(:,2,2) = dy2dx2; % = dftdx1-dfbdx1/2
                % y2 = f_b(x1) + ( (x2+1)/2 ) * (ft(x1) - fb(x1))
                J(:,2,1) = dfbdx1 + (dftdx1-dfbdx1).*(x2+1)/2;
            end

            if(nargout >= 4)
                % Hessian of y1 map
                dH1        = zeros(n,2,2);
            end

            if(nargout >= 4)
                % Hessian of y2 map
                % y2 = bottom(x1) + (x2+1).*(top(x1)-bottom(x1))/2;
                dH2        = zeros(n,2,2);
                dH2(:,1,1) = dfbdxx1 + (dftdxx1-dfbdxx1).*(x2+1)/2;
                dH2(:,1,2) = (dftdx1-dfbdx1)/2;
                dH2(:,2,1) = (dftdx1-dfbdx1)/2;
            end

        end

        function CompSpace(~,~,~)
            % No inverse map yet
            exc = MException('FluctuatingChannel:CompSpace','not yet implemented');
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

            % shorthand parameters
            y1Min = 0; 
            L1 = this.L1; %#ok<*NASGU>
            y1Max = this.L1;

            % create top and bottom on a fine grid for plotting
            xPlot = (-1:1/100:1)';
            y1Plot = LinearMap(xPlot,y1Min,y1Max);
            top = this.ft(xPlot,this.ftParams);
            bottom = this.ft(xPlot,this.fbParams);
          
            % shorthand parameters
            y2Min = min(bottom);
            y2Max = max(top);
            L2 = y2Max - y2Min;

            % set sensible axis limits
            xLim = [y1Min - 0.1*L1, y1Max + 0.1*L1];
            yLim = [y2Min - 0.1*L2, y2Max + 0.1*L2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;

            % plot the boundary
            hold on
            plot([y1Min,y1Min],[bottom(1),top(1)],'-r','DisplayName', 'Boundary');
            plot([y1Max,y1Max],[bottom(end),top(end)],'-r','HandleVisibility','off');
            plot(y1Plot,bottom,'-r','HandleVisibility','off');
            plot(y1Plot,top,'-r','HandleVisibility','off');

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

            geom.N = [20;20];
            geom.L1 = 3;
            geom.ft='ChannelSin'; 
            geom.fb = 'ChannelSin';
            ftParams.shift = 1; 
            fbParams.shift = 0; 
            ftParams.omega = 2*pi/geom.L1; 
            fbParams.omega = 2*pi/geom.L1;
            ftParams.A = 0.2; 
            fbParams.A = -0.2;
            ftParams.phi = 0; 
            fbParams.phi = 0;
            geom.ftParams = ftParams; 
            geom.fbParams = fbParams;

        end

    end


end


