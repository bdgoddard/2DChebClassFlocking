classdef HalfSpaceMinusDisk < ComposedShape
    % An upper half space without a disc of radius R at Origin.
    % This is a ComposedShape.
    % Note that the wall, given by y2Wall, must be below the Origin of the
    % disc.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R - radius of the disc to be removed
    %  y2Wall - location of the wall
    %  L1 - mapping parameter for radial coordinate; half the points lie in
    %       [R,R+L1]
    %  (Origin = [0,0]) - origin of the disc


    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = HalfSpaceMinusDisk(Geometry)

            % create the half space without half a disc
            this.SubShape{1} = HalfSpaceMinusHalfDisk(Geometry);

            % deal with the remaining part
            if((Geometry.Origin(2)-Geometry.y2Wall) > Geometry.R)
                % if the whole disc is in the half space
                this.SubShape{2} = StripMinusDisk(Geometry);

            elseif(((Geometry.Origin(2)-Geometry.y2Wall) < Geometry.R) ...
                    && ((Geometry.Origin(2)-Geometry.y2Wall) > 0))

                % if only part of the disc is in the half space
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                Geometry.TopBottom = 'Bottom';
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);

                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);

            elseif(((Geometry.Origin(2)-Geometry.y2Wall) < 0))
                % if the Origin is below the wall
                error('Not yet implemented');
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
            L1 = this.SubShape{1}.L1;
            R  = this.SubShape{1}.R;
            y2Wall = this.SubShape{2}.y2Wall;
            O1 = this.SubShape{1}.Origin(1);
            O2 = this.SubShape{1}.Origin(2);
            deep = (O2 - y2Wall) < R;
            if(deep)
                th2 = this.SubShape{2}.th2;
            end
            L1Strip = abs(O2 - y2Wall);


            % set sensible axis limits
            xLim = [-1.1 * 3*L1 + O1, 1.1 * 3*L1 + O1];
            yLim = [y2Wall - 0.1*3*L1Strip, 1.1 * 3*L1 + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grid
            this.PlotGrid;
            hold on

            % plot indicator line for L1+R
            tCircle = linspace(0, 2*pi, 100);
            plot((L1+R)*cos(tCircle) + O1, (L1+R)*sin(tCircle) + O2, ...
                '--m', 'DisplayName','L1+R');

            % plot indicator line for domain split
            plot(xLim,[O2,O2],':k', 'DisplayName','Domain split');


            % plot circle
            tCircle = linspace(0, 2*pi, 100);
            plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                '--k', 'DisplayName','Wall and Disc');
            p1 = [O1, O2];
            theta = 3*pi/4;
            p2 = [O1 + R * cos(theta), O2 + R * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)], '-b','HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * 3*L1, p1(2) + dp(2)/2 + 0.02 * 3*L1, ['R = ' num2str(R)], ...
                'FontSize',20,'Color','b')

            % plot wall
            plot(xLim,[y2Wall,y2Wall],'--k', 'HandleVisibility','off');

            if(deep)
                % plot boundary lines
                plot([O1 - R*cos(th2),xLim(2)],[y2Wall,y2Wall],'-r', 'DisplayName', 'Boundary');
                plot([xLim(1),O1 + R*cos(th2)],[y2Wall,y2Wall],'-r', 'HandleVisibility','off');
                % plot boundary arc
                tCircle = linspace(th2, pi-th2, 100);
                plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                    '-r', 'HandleVisibility','off');
            else
                % plot boundary line
                plot(xLim,[y2Wall,y2Wall],'-r', 'DisplayName', 'Border');
                % plot circle
                tCircle = linspace(0, 2*pi, 100);
                plot(R*cos(tCircle) + O1, R*sin(tCircle) + O2, ...
                    '-r', 'HandleVisibility','off');
            end

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',5,'MarkerFaceColor','k','DisplayName','Origin')

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

            geom.N = [20;12];
            geom.R = 1;
            geom.L1 = 1;
            % geom.y2Wall = -0.6; % option 1
            geom.y2Wall = -1.2; % option 2
            geom.Origin = [0,0];

        end

    end

end

