classdef AnnulusCut < ComposedShape
    % An annulus with innner radius R_in and outer radius R_in cut off by a
    % wall that is h below the Origin.  This is a ComposedShape with
    % SubShape structure.
    %
    % Inputs:
    %  N - [N1,N2] number of points in the two coordinates
    %  R_out - outer radius of the annulus
    %  R_in  - inner radius of the annulus
    %  h - perpendicular distance from the Origin to the wall
    %  (Origin = [0,0]) - origin of the annulus

    methods

        %------------------------------------------------------------------
        % Constructor
        %------------------------------------------------------------------

        function this = AnnulusCut(Geometry)

            % determine some useful angles
            % these are where the outer cicle meets the wall
            th             = acos(Geometry.h/Geometry.R_out);
            Geometry.th1   = 3/2*pi + th;
            Geometry.th2   = 3/2*pi - th;

            % a Wedge for the part that is not cut off
            this.SubShape{1}  = Wedge(Geometry);

            if((Geometry.h > 0) && (Geometry.h < Geometry.R_in))
                % cut off across both circles
                % results in two separate wedges
                Geometry.leftRight = 'left';
                this.SubShape{2} = WedgeCutSide(Geometry);

                Geometry.leftRight = 'right';
                this.SubShape{3} = WedgeCutSide(Geometry);

            elseif((Geometry.h > 0) && (Geometry.h < Geometry.R_out))
                % cut off only through the outer circle
                this.SubShape{2} = WedgeCut(Geometry);

            elseif((Geometry.h >= Geometry.R_out))
                % not cut off at all
                error('AnnulusCut: Define full annulus')
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
            R_in = this.SubShape{2}.R_in;
            R_out = this.SubShape{2}.R_out;
            h = this.SubShape{2}.h;
            th1 = this.SubShape{2}.th1;
            th2 = this.SubShape{2}.th2;
            O1 = this.SubShape{2}.Origin(1);
            O2 = this.SubShape{2}.Origin(2);

            % set sensible axis limits
            xLim = [-1.1 * R_out + O1, 1.1 * R_out + O1];
            yLim = [-1.1 * R_out + O2, 1.1 * R_out + O2];

            % set up figure
            hf = figure('Position',[100,100,1200,800]);
            ha = axes;

            % plot the grd
            this.PlotGrid;
            hold on

            % plot the origin
            plot(O1,O2,'ok','MarkerSize',15,'MarkerFaceColor','k','DisplayName','Origin')

            % plot inner and outer circles
            % R_in
            if(R_in > 0)
                tCircle = linspace(0, 2*pi, 100);
                plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                    '--k', 'HandleVisibility','off');
                p1 = [O1, O2];
                theta = 3*pi/4;
                p2 = [O1 + R_in * cos(theta), O2 + R_in * sin(theta)];
                dp = p2 - p1;
                plot([p1(1),p2(1)],[p1(2),p2(2)], '-b', 'HandleVisibility','off');
                text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2 + 0.02 * R_out, 'R in', ...
                    'FontSize',20,'Color','b')
            end

            % R_out
            tCircle = linspace(0, 2*pi, 100);
            plot(R_out*cos(tCircle) + O1, R_out*sin(tCircle) + O2, ...
                    '--k', 'DisplayName','Circle and wall');
            p1 = [O1, O2];
            theta = pi/4;
            p2 = [O1 + R_out * cos(theta), O2 + R_out * sin(theta)];
            dp = p2 - p1;
            plot([p1(1),p2(1)],[p1(2),p2(2)], '-b', 'HandleVisibility','off');
            text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2 - 0.02 * R_out, 'R out', ...
                'FontSize',20,'Color','b')

            % plot wall 
            plot(xLim,[O2 - h, O2 - h],'--k', 'HandleVisibility','off');


            % plot indicators for the two angles
            tth1 = linspace(0, mod(th1,2*pi), 100);
            tth2 = linspace(0, mod(th2,2*pi), 100);
            R1 = R_in + (R_out-R_in)/4;
            R2 = R_in + (R_out-R_in)/3;
            plot(R1 * cos(tth1) + O1, R1 * sin(tth1) + O2, ...
                '-m', 'DisplayName', 'th1, th2');
            plot(R2 * cos(tth2) + O1, R2 * sin(tth2) + O2, ...
                '-m', 'HandleVisibility','off')
            % plot lines at theta = 0 and theta = th2
            plot([O1,xLim(2)],[O2,O2],':m','HandleVisibility','off')
            plot([R_in*cos(th2) + O1, R_out*cos(th2) + O1],[R_in*sin(th2) + O2, R_out*sin(th2) + O2], ...
                        ':m', 'HandleVisibility','off');

            % plot indicator for h
            p1 = [O1, O2];
            p2 = [O1, O2 - h];
            dp = p2 - p1;
            plot([p1(1), p2(1)], [p1(2), p2(2)], '-b', 'DisplayName','h');
            text((p1(1) + dp(1)/2) + 0.02 * R_out, p1(2) + dp(2)/2, 'h', ...
                'FontSize',20,'Color','b')

            % plot the boundary
            % outer circle
            tCircle = linspace(th1, pi - th1, 100);
            plot(R_out*cos(tCircle) + O1, R_out*sin(tCircle) + O2, ...
                    '-r', 'DisplayName','Boundary');

            if( h < R_in) % double cut

                % plot the inner circle
                if(R_in > 0)
                    tCircle = linspace(th2, pi - th2, 100);
                    plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                        '-r', 'HandleVisibility','off');
                end

                % bottom lines
                plot([R_out*cos(th1) + O1, R_in*cos(th2) + O1],[O2 - h, O2 - h], ...
                        '-r', 'HandleVisibility','off');
                plot([R_out*cos(pi-th1) + O1, R_in*cos(pi-th2) + O1],[O2 - h, O2 - h], ...
                        '-r', 'HandleVisibility','off');

                % plot domain decompositions
                plot([R_in*cos(th1) + O1, R_out*cos(th1) + O1],[R_in*sin(th1) + O2, R_out*sin(th1) + O2], ...
                        ':k', 'DisplayName','Domain split');
                plot([R_in*cos(pi-th1) + O1, R_out*cos(pi-th1) + O1],[R_in*sin(pi-th1) + O2, R_out*sin(pi-th1) + O2], ...
                        ':k', 'HandleVisibility','off');


            else

                % plot the inner circle
                if(R_in > 0)
                    tCircle = linspace(0, 2*pi, 100);
                    plot(R_in*cos(tCircle) + O1, R_in*sin(tCircle) + O2, ...
                        '-r', 'HandleVisibility','off');
                end

                % bottom line
                plot([R_out*cos(th1) + O1, R_out*cos(th2) + O1],[R_out*sin(th2) + O2, R_out*sin(th2) + O2], ...
                        '-r', 'HandleVisibility','off');

                % plot domain decomposition
                plot([R_in*cos(th1) + O1, R_out*cos(th1) + O1],[R_in*sin(th1) + O2, R_out*sin(th1) + O2], ...
                        ':k', 'DisplayName','Domain split');
                plot([R_in*cos(th2) + O1, R_out*cos(th2) + O1],[R_in*sin(th2) + O2, R_out*sin(th2) + O2], ...
                        ':k', 'HandleVisibility','off');


            end

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
            geom.R_in = 2;
            geom.R_out = 5;
            % geom.h = 2.5; % case 1
            geom.h = 1.8; % case 2
            geom.Origin = [-0.2,0.5];

        end

    end



end

