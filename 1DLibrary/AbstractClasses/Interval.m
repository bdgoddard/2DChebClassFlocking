classdef (Abstract) Interval < handle
    
    properties (Access = public)
        Pts
        Diff
        Int
        Ind
        Interp
        Conv
        
        N
        polar = 'undefined';
    end
    
    methods (Abstract=true, Access=public)
        [y,J,dH,ddH] = PhysSpace(this,x);
        [x]      = CompSpace(this,y);
        
        Ind    = ComputeIndices(this);
        Diff   = ComputeDifferentiationMatrix(this);
        Interp = ComputeInterpolationMatrix(this,interp,saveBool);           
        Int    = ComputeIntegrationVector(this);
        M_conv = ComputeConvolutionMatrix(this,f,saveBool);
    end
    
    methods
        function this = Interval(N)
            this.N=N;
        end
    end
    
    methods (Access = public)
        
        function d  = GetDistance(this,pts_y1,pts_y2)                         
            d       = abs(pts_y1 - pts_y2);            
        end        
        function [Pts,Diff,Int,Ind,Interp,Conv] = ComputeAll(this,PlotRange,f)
            Ind  = ComputeIndices(this);
            Diff = ComputeDifferentiationMatrix(this);
            Int  = ComputeIntegrationVector(this);
                        
            Pts    = this.Pts;            
            
            if(nargin >= 2)
                Interp = InterpolationPlot(this,PlotRange,true);  
            end
            
            if(nargin >= 3)
                    Conv = ComputeConvolutionMatrix(this,f,true);
            end
        end                
        function [IP,yPlot] = InterpolationPlot(this,PlotRange,saveBool)
            
            xMin = CompSpace(this,PlotRange.yMin);
            xMax = CompSpace(this,PlotRange.yMax);
            xPlot = GetArray(xMin,xMax,PlotRange.N);

            yPlot = PhysSpace(this,xPlot);
            
            if(nargin == 3)            
                IP = ComputeInterpolationMatrix(this,xPlot,saveBool);           
            else
                IP = ComputeInterpolationMatrix(this,xPlot,false);           
            end
        end  
        
        function IP = ComputeInterpolationMatrixPhys(this,y,saveBool)
                        
            x = CompSpace(this,y);
            
            if(nargin == 3)            
                IP = ComputeInterpolationMatrix(this,x,saveBool);           
            else
                IP = ComputeInterpolationMatrix(this,x,false);           
            end
        end  
        
        function InterPolUniform = InterpolationPlotUniform(this,dy)   

            if(~isfield(this.Pts,'Uniform'))
                if(nargin<2)
                    this.ComputeUniformPts;
                else
                    this.ComputeUniformPts(dy);
                end
            end
                        
            this.Interp.InterPolUniform = this.InterpolationMatrix_Pointwise(this.Pts.Uniform.y);
            
            InterPolUniform = this.Interp.InterPolUniform;
            
        end

        function UniformPts = ComputeUniformPts(this,dy)
            
            % find bounds            
            yMin = min(this.Pts.y);
            yMax = max(this.Pts.y);
            
            % find uniform step
            if(nargin<2)
                L = yMax-yMin;
                NPts = 10;
                dy = L/NPts;
            end
            
            % find bounds for uniform grid
            yMin0 = floor(yMin/dy)*dy;
            yMax0 = ceil(yMax/dy)*dy;
            
            % construct uniform grid
            y = (yMin0:dy:yMax0)';
            x = this.CompSpace(y);
            
            % remove those outside the shape (i.e. which get mapped outside
            % [-1,1] x [-1,1] in the computational space
            outOfBounds = (x > 1) | (x < -1);
            y(outOfBounds) = [];
            x(outOfBounds) = [];
            
            this.Pts.Uniform = struct('y',y,'x',x);
            
            UniformPts = this.Pts.Uniform;
            
        end

        
        
        function InitializationPts(this)
            this.Pts.N = this.N;
            this.Pts.y = PhysSpace(this,this.Pts.x);
        end           
        
        function h = plot(this,V,options)
            
            nSpecies = size(V,2);
            
            if(nSpecies==1)
                nCol = 1;
            else
                nCol = 2;
            end            
            nRows = ceil(nSpecies/nCol);
            
            if((nargin > 2) && strcmp(options,'plain'))
                options = struct('plain',true,...
                                 'linecolor','k','linewidth',1.5);
            elseif(nargin<3)
                options = struct;
            end
            
            if(isfield(options,'ha'))
                ha = options.ha;
            else
                ha = gca;
            end
                        
            for iSpecies = 1:nSpecies
                if(nSpecies>1)
                    subplot(nRows,nCol,iSpecies)
                end
                
                if( (size(V,1) == length(this.Interp.pts))  ...
                        && (length(this.Interp.pts) ~= length(this.Pts.y)) )
                    VI = V(:,iSpecies);
                else
                    VI = this.Interp.InterPol*V(:,iSpecies);
                end

                h = plot(ha,this.Interp.pts,VI);
                
                if(nargin>2)                    
                    if(isfield(options,'linecolor'))
                        set(h,'Color',options.linecolor);
                    end
                    if(isfield(options,'linestyle'))
                        set(h,'LineStyle',options.linestyle);
                    end
                    if(isfield(options,'linewidth'))
                        set(h,'LineWidth',options.linewidth);
                    end
                    if(isfield(options,'marker'))
                        set(h,'Marker',options.marker);
                        set(h,'MarkerSize',15);
                    end                    
                    if(isfield(options,'xLim'))
                        xlim(options.xLim);
                    end                    
                    if(isfield(options,'yLim'))
                        ylim(options.yLim);
                    end
                    if(isfield(options,'displayname'))
                        if(iscell(options.displayname))
                            options.displayname = options.displayname{1};
                        end
                        set(h,'DisplayName',options.displayname);
                    end
                    if(isfield(options,'handlevisibility'))
                        set(h,'HandleVisibility',options.handlevisibility);
                    end

                end                
                
                if((nargin <= 2) || ~isfield(options,'plain') || ~options.plain)
                    hold on;
                    plot(this.Pts.y,V,'o','MarkerEdgeColor','k','MarkerFaceColor','g');
                end
                
                if(nargin>2)
                    if(isfield(options,'xlabel'))
                        xlabel(options.xlabel,'Interpreter','Latex','fontsize',25);
                    else
                        xlabel('$y$','Interpreter','Latex','fontsize',25);
                    end
                    if(isfield(options,'ylabel'))
                        ylabel(options.ylabel,'Interpreter','Latex','fontsize',25);
                    else
                        ylabel('$\rho$','Interpreter','Latex','fontsize',25);
                    end
                end
                
                xlim([min(this.Interp.pts) max(this.Interp.pts)]);
                
                if(nSpecies > 1)
                    title(['Species ' num2str(iSpecies)]);
                end
                
                set(ha,'fontsize',20);
                set(ha,'linewidth',1.5);    
                
            end
        end        
        function PlotGrid(this)
            hold on
            %yRange = get(gca,'YLim');
            %h = plot(this.Pts.y, yRange(1),'og');
            h = plot(this.Pts.y, zeros(size(this.Pts.y)),'og');
            set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold off
        end        
        function IP   = InterpolationMatrix_Pointwise(this,yP)
            IP = zeros(length(yP),this.N);            
            for i =1:length(yP)
                x       = CompSpace(this,yP(i));
                h       = ComputeInterpolationMatrix(this,x,false);
                IP(i,:) = h.InterPol;
            end
        end        
        function IP = SubShapePts(this,a_shapePts)                   
            IP = InterpolationMatrix_Pointwise(this,a_shapePts.y);                           
        end      
        function IntM = IntegrationY(this)
            
            y    = this.Pts.y;
            IntM = zeros(length(y));
            vh   = zeros(1,length(y));
            for i = 2:length(y)
                %Matrix for integral int(v(y),y=y1..Y)
                h           = y(i) - y(i-1);
                if(h == inf)   
                    h = 0; %assume value to be integrated is zero, if h = inf
                end
                vh([i-1,i]) = vh([i-1,i]) + h/2;
                
                IntM(i,:) = vh;
            end
            

        end
    end
    
    
end