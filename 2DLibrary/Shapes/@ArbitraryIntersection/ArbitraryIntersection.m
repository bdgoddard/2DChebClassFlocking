classdef ArbitraryIntersection < matlab.mixin.Copyable
    % acts as a handle class but allows copying via
    % copy(ArbitraryIntersection), otherwise it's the same pointer and
    % shifting, etc overwrites the original version
    
    properties
        Pts
        Int
        Area
        Origin = [0;0];
    end
    
    methods
        function this = ArbitraryIntersection(Shapes)
            if(nargin == 0)
                nShapes = 0;
            else
                nShapes = length(Shapes);
            end
            
            y1_kv = [];
            y2_kv = [];
            this.Int = [];
            this.Area = 0;
            
            for iShape = 1:nShapes
                shapePts = Shapes(iShape).Shape.GetCartPts;
                y1_kv = [y1_kv ; shapePts.y1_kv];
                y2_kv = [y2_kv ; shapePts.y2_kv];
                [int,area] = Shapes(iShape).Shape.ComputeIntegrationVector;
                this.Int = [this.Int , int];
                this.Area = this.Area + area;
            end
            
            this.Pts.y1_kv = y1_kv;
            this.Pts.y2_kv = y2_kv;
            
            if(nargin>0)
                this.Origin = Shapes(1).Shape.Origin;
            end
        end
        
        
        function ptsCart = GetCartPts(this)
            ptsCart = this.Pts;
        end
        
        
        function [int,area] = ComputeIntegrationVector(this)
            int = this.Int;
            area = this.Area;
        end
        
        function ShiftPts(this,shiftVec)
            this.Pts.y1_kv = this.Pts.y1_kv + shiftVec(1);
            this.Pts.y2_kv = this.Pts.y2_kv + shiftVec(2);
        end
        
        function RotatePts(this,theta)
            y1_kv = this.Pts.y1_kv;
            y2_kv = this.Pts.y2_kv;
            this.Pts.y1_kv = y1_kv*cos(theta) - y2_kv*sin(theta); % can have issues when inf - inf = nan
            this.Pts.y2_kv = y1_kv*sin(theta) + y2_kv*cos(theta);
            this.Pts.y1_kv(isnan(this.Pts.y1_kv)) = inf; % need to think about this some more
            this.Pts.y2_kv(isnan(this.Pts.y2_kv)) = inf;
        end
        
        function ReflectPtsXAxis(this)
            this.Pts.y2_kv = - this.Pts.y2_kv;
        end
        
        function ReflectPtsYAxis(this)
            this.Pts.y1_kv = - this.Pts.y1_kv;
        end
        
        
        function h = PlotGrid(this,optsPlot)
            
            if(nargin<2)
                optsPlot = struct;
            end
            
            facecolor = 'g';
            edgecolor = 'k';
            marker = 'o';
            
            if(isfield(optsPlot,'facecolor'))
                facecolor = optsPlot.facecolor;
            end
            if(isfield(optsPlot,'edgecolor'))
                edgecolor = optsPlot.edgecolor;
            end
            if(isfield(optsPlot,'marker'))
                marker = optsPlot.marker;
            end
            
            ptsCart = GetCartPts(this);
            h = scatter(ptsCart.y1_kv,ptsCart.y2_kv);
            set(h, 'Marker', marker, ...
                'MarkerEdgeColor',edgecolor,...
                'MarkerFaceColor',facecolor );
            if(isfield(optsPlot,'displayname'))
                set(h,'DisplayName',optsPlot.displayname);
            else
                set(h,'HandleVisibility','off');
            end
            
        end
        
        
    end
    
    
    
end