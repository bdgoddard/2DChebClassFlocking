classdef ArbitraryIntersection < matlab.mixin.Copyable
    % Acts as a handle class but allows copying via
    % copy(ArbitraryIntersection), otherwise it's the same pointer and
    % shifting, etc overwrites the original version
    %
    % Designed to take a set of Shapes and combine the points (Pts),
    % integration weights (Int), and areas (Area).
    % Not suitable for general computations, but can do integration.
    %
    % Inputs:
    %  Shapes - a structure of Shapes, which must in turn contain:
    %    Origin - the origin
    %    GetCartPts - a method to return the Cartesian points as a structure
    %                 containing y1_kv and y2_kv
    %    ComputeIntegrationVector - a method returning the integration weights
    %                 as a row vector and the area as a scalar
    %

    properties
        Pts % points
        Int % integration vector
        Area % exact area (requires exact area from Shapes)
        Mask % store where each set of points came from
        Origin = [0;0];  % origin
    end

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = ArbitraryIntersection(Shapes)

            % can construct an empty ArbitraryIntersection
            if(nargin == 0)
                nShapes = 0;
            else
                nShapes = length(Shapes);
            end

            % set up initial values of points, integration vector, and area
            y1_kv = [];
            y2_kv = [];
            this.Int = [];
            this.Area = 0;
            this.Mask = [];

            % loop through the Shapes and concatenate the points,
            % integration vectors, and areas
            for iShape = 1:nShapes

                shapePts = Shapes(iShape).Shape.GetCartPts;
                y1_kv = [y1_kv ; shapePts.y1_kv]; %#ok<*AGROW>
                y2_kv = [y2_kv ; shapePts.y2_kv];

                [int,area] = Shapes(iShape).Shape.ComputeIntegrationVector;

                % Int mimics the form in Shape
                this.Int = [this.Int , int];

                this.Area = this.Area + area;

                % make a mask to store where each point came from
                this.Mask = [this.Mask; iShape*ones(size(shapePts.y1_kv))];
            end

            % set the Pts structure to mimic Shape
            this.Pts.y1_kv = y1_kv;
            this.Pts.y2_kv = y2_kv;

            if(nargin>0)
                this.Origin = Shapes(1).Shape.Origin;
            end
        end

        %------------------------------------------------------------------
        % Shape-like functions
        %------------------------------------------------------------------

        function ptsCart = GetCartPts(this)
            % mimic Shape functionality
            ptsCart = this.Pts;
        end


        function [int,area] = ComputeIntegrationVector(this)
            % mimic Shape functionality
            int = this.Int;
            area = this.Area;
        end

        function h = PlotGrid(this,optsPlot)
            % mimic Shape's plotting of grid

            if(isempty(this.Pts.y1_kv))
                h = [];
                return
            end

            % set some options
            if(nargin<2)
                optsPlot = struct;
            end

            % plot by subshape
            if(isfield(optsPlot,'subshapes'))
                subshapes = optsPlot.subshapes;
                nSubshapes = this.Mask(end);
            else
                subshapes = false;
            end

            % default plotting to match Shape
            facecolor = 'g';
            edgecolor = 'k';
            marker = 'o';
            linewidth = 0.5;
            
            if(isfield(optsPlot,'facecolor'))
                facecolor = optsPlot.facecolor;
            end

            if(isfield(optsPlot,'edgecolor'))
                edgecolor = optsPlot.edgecolor;
            end

            if(isfield(optsPlot,'marker'))
                marker = optsPlot.marker;
            end

            if(isfield(optsPlot,'linewidth'))
                linewidth = optsPlot.linewidth;
            end


            % plot the points
            ptsCart = GetCartPts(this);
            if(subshapes)
                h = zeros(nSubshapes,1);
                for iSubshape = 1:nSubshapes
                    iMask = this.Mask == iSubshape;

                    if(iscell(facecolor) && length(facecolor) >= nSubshapes)
                        iFacecolor = facecolor{iSubshape};
                    elseif(~iscell(facecolor) && size(facecolor,1) >= nSubshapes)
                        iFacecolor = facecolor(iSubshape,:);
                    else
                        iFacecolor = facecolor(1,:);
                    end

                    if(iscell(edgecolor) && length(edgecolor) >= nSubshapes)
                        iEdgecolor = edgecolor{iSubshape};
                    elseif(~iscell(edgecolor) && size(edgecolor,1) >= nSubshapes)
                        iEdgecolor = edgecolor(iSubshape,:);
                    else
                        iEdgecolor = edgecolor(1,:);
                    end

                    if(length(marker) >= nSubshapes)
                        iMarker = marker{iSubshape};
                    else
                        iMarker = marker{1};
                    end

                    if(length(linewidth) >= nSubshapes)
                        iLinewidth = linewidth(iSubshape);
                    else
                        iLinewidth = linewidth(1);
                    end

                    h(iSubshape) = scatter(ptsCart.y1_kv(iMask),ptsCart.y2_kv(iMask), ...
                                            'Marker', iMarker, ...
                                            'MarkerFaceColor',iFacecolor, ...
                                            'MarkerEdgeColor',iEdgecolor, ...
                                            'LineWidth',iLinewidth);
                    
                    if(isfield(optsPlot,'size'))
                        set(h(iSubshape),'SizeData',optsPlot.size);
                    end

                    hold on

                end
            else
                h = scatter(ptsCart.y1_kv,ptsCart.y2_kv);

                % make the plot nice
                set(h, 'Marker', marker, ...
                    'MarkerEdgeColor',edgecolor,...
                    'MarkerFaceColor',facecolor, ...
                    'LineWidth',linewidth);
                if(isfield(optsPlot,'size'))
                    set(h,'SizeData',optsPlot.size);
                end

                % allow for using a legend
                if(isfield(optsPlot,'displayname'))
                    set(h,'DisplayName',optsPlot.displayname);
                else
                    set(h,'HandleVisibility','off');
                end


            end
            
            

        end

        %------------------------------------------------------------------
        % Geometrical functions
        %------------------------------------------------------------------

        function ShiftPts(this,shiftVec)
            % shift points by shiftVec
            this.Pts.y1_kv = this.Pts.y1_kv + shiftVec(1);
            this.Pts.y2_kv = this.Pts.y2_kv + shiftVec(2);
        end

        function RotatePts(this,theta)
            % rotate by theta around (0,0)

            % shorthand the points
            y1_kv = this.Pts.y1_kv;
            y2_kv = this.Pts.y2_kv;

            % do the rotation
            this.Pts.y1_kv = y1_kv*cos(theta) - y2_kv*sin(theta); % can have issues when inf - inf = nan
            this.Pts.y2_kv = y1_kv*sin(theta) + y2_kv*cos(theta);
            
            % (attempt to) deal with rotation when some points are inf
            this.Pts.y1_kv(isnan(this.Pts.y1_kv)) = inf; % need to think about this some more
            this.Pts.y2_kv(isnan(this.Pts.y2_kv)) = inf;
        end

        function ReflectPtsXAxis(this)
            % reflect in the x-axis
            this.Pts.y2_kv = - this.Pts.y2_kv;
        end

        function ReflectPtsYAxis(this)
            % reflect in the y-axis
            this.Pts.y1_kv = - this.Pts.y1_kv;
        end

    end



end