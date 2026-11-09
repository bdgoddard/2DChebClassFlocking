classdef (Abstract) Shape < handle

    properties (Access = public)
        Pts
        Diff
        Int
        Ind
        Interp
        InterpFlux = []
        Conv
        Origin = [0,0]

        BoundaryPaths %vector with [Path_right; PathUpper ; PathLeft ; PathBottom]

        N1
        N2
        M
        polar = 'undefined';

        periodic = [false;false];

        invertibleMap = true; % true if the map is invertible.  If not, we can't differentiate.
    end

    methods (Abstract = true,Access = public)
        [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(x1,x2);
        [x1,x2]                 = CompSpace(y1,y2);

        Ind    = ComputeIndices(this);
        Diff   = ComputeDifferentiationMatrix(this);
        Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool);
        Int    = ComputeIntegrationVector(this);
        M_conv = ComputeConvolutionMatrix(this,f,saveBool);
    end

    methods
        function this = Shape(N1,N2)
            this.N1 = N1;
            this.N2 = N2;
        end
    end

    methods (Access = public)

        function ptsCart = GetCartPts(this,pts_y1,pts_y2)

            if(nargin == 1)
                pts_y1 = this.Pts.y1_kv;
                pts_y2 = this.Pts.y2_kv;
            end

            if(strcmp(this.polar,'polar'))
                [ptsCart.y1_kv,ptsCart.y2_kv] = pol2cart(pts_y2,pts_y1);
                mark0 = (((pts_y2==2*pi) |(pts_y2==0) | (pts_y2==pi)) & (pts_y1==Inf));
                ptsCart.y2_kv(mark0) = 0;
            elseif(strcmp(this.polar,'sphSurf'))
                th  = pts_y1;  phi = pts_y2;
                th0 = this.theta0;
                y1_kv_temp = this.R*sin(th).*cos(phi);
                y2_kv_temp = this.R*cos(th);
                ptsCart.y1_kv = cos(th0)*y1_kv_temp - sin(th0)*y2_kv_temp;
                ptsCart.y2_kv = sin(th0)*y1_kv_temp + cos(th0)*y2_kv_temp;
            elseif(strcmp(this.polar,'cart'))
                ptsCart.y1_kv = pts_y1;
                ptsCart.y2_kv = pts_y2;
            else
                exc = MException('Shape:GetCartPts','select {polar,sphSurf,cart}');
                throw(exc);
            end

            ptsCart.y1_kv = ptsCart.y1_kv + this.Origin(1);
            ptsCart.y2_kv = ptsCart.y2_kv + this.Origin(2);

        end

        function pts = GetInvCartPts(this,ptsCart_y1,ptsCart_y2)

            if(strcmp(this.polar,'polar'))
                [pts.y2_kv,pts.y1_kv] = cart2pol(ptsCart_y1,ptsCart_y2);
            elseif(strcmp(this.polar,'sphSurf'))
                %th  = pts_y1;  phi = pts_y2;
                %ptsCart.y1_kv = this.R*sin(th).*cos(phi);
                %ptsCart.y2_kv = this.R*cos(th);
                exc = MException('Shape:GetInvCartPts','sphSurf not yet implemented');
                throw(exc);
            elseif(strcmp(this.polar,'cart'))
                pts.y1_kv = ptsCart_y1;
                pts.y2_kv = ptsCart_y2;
            else
                exc = MException('Shape:GetInvCartPts','select {polar,sphSurf,cart}');
                throw(exc);
            end

        end

        function vecPtsCart = GetCartVec(this,vecPts)
            % Get Cartesian representation of a list of r and theta values
            % treated as determining a vector x, y

            if(nargin == 1)
                tempPts_y1 = this.Ind.normal(:,1:end/2);
                tempPts_y2 = this.Ind.normal(:,end/2+1:end);
            else
                rowVec = ( size(vecPts,2) == 2*this.M );
                if(rowVec)  % row vectors
                    tempPts_y1 = vecPts(:,1:end/2);
                    tempPts_y2 = vecPts(:,end/2+1:end);
                else  % column vectors
                    tempPts_y1 = vecPts(1:end/2,:);
                    tempPts_y2 = vecPts(end/2+1:end,:);
                end
            end

            if(strcmp(this.polar,'polar'))
                [tempPtsCart.y1_kv,tempPtsCart.y2_kv] = pol2cart(tempPts_y2,tempPts_y1);
                mark0 = (((tempPts_y2==2*pi) |(tempPts_y2==0) | (tempPts_y2==pi)) & (tempPts_y1==Inf));
                tempPtsCart.y2_kv(mark0) = 0;
            elseif(strcmp(this.polar,'sphSurf'))
                th  = tempPts_y1;  phi = tempPts_y2;
                tempPtsCart.y1_kv = this.R*sin(th).*cos(phi);
                tempPtsCart.y2_kv = this.R*cos(th);
            elseif(strcmp(this.polar,'cart'))
                tempPtsCart.y1_kv = tempPts_y1;
                tempPtsCart.y2_kv = tempPts_y2;
            else
                exc = MException('Shape:GetCartVec','select {polar,sphSurf,cart}');
                throw(exc);
            end

            % tempPtsCart.y1_kv = tempPtsCart.y1_kv + this.Origin(1);
            % tempPtsCart.y2_kv = tempPtsCart.y2_kv + this.Origin(2);

            if(rowVec)
                vecPtsCart = [tempPtsCart.y1_kv, tempPtsCart.y2_kv];
            else
                vecPtsCart = [tempPtsCart.y1_kv; tempPtsCart.y2_kv];
            end

        end

        function d  = GetDistance(this,pts_y1,pts_y2)
            if(nargin == 2)
                pts_y2 = pts_y1.y2_kv;
                pts_y1 = pts_y1.y1_kv;
            end

            ptsCart = GetCartPts(this,pts_y1,pts_y2);
            d       = sqrt(ptsCart.y1_kv.^2 + ptsCart.y2_kv.^2);
        end

        function [Pts,Diff,Int,Ind,Interp,Conv] = ComputeAll(this,PlotArea,f,shapeParam)

            Ind  = ComputeIndices(this);
            Diff = ComputeDifferentiationMatrix(this);
            Int  = ComputeIntegrationVector(this);

            Pts    = this.Pts;

            if(nargin >= 2)
                Interp = InterpolationPlot(this,PlotArea,true);
            end
            if(nargin >= 3)
                if(nargin ==4)
                    Conv = ComputeConvolutionMatrix(this,f,shapeParam);
                else
                    Conv = ComputeConvolutionMatrix(this,f,true);
                end
            end
        end

        function this = InitializationPts(this)
            this.Pts.x1_kv  = kron(this.Pts.x1,ones(this.N2,1));
            this.Pts.x2_kv  = kron(ones(this.N1,1),this.Pts.x2);

            [this.Pts.y1_kv,this.Pts.y2_kv] = PhysSpace(this,this.Pts.x1_kv,...
                this.Pts.x2_kv);
            this.Pts.N1 = this.N1;
            this.Pts.N2 = this.N2;
            this.M      = length(this.Pts.y1_kv);
        end

        function UniformPts = ComputeUniformPts(this,dy)

            tempPts = this.Pts;

            if(strcmp(this.polar,'polar'))
                ptsCart = this.GetCartPts(this.Pts.y1_kv,this.Pts.y2_kv);
                this.Pts.y1_kv = ptsCart.y1_kv;
                this.Pts.y2_kv = ptsCart.y2_kv;
            end
            % find bounds
            y1Min = min(this.Pts.y1_kv);
            y1Max = max(this.Pts.y1_kv);
            y2Min = min(this.Pts.y2_kv);
            y2Max = max(this.Pts.y2_kv);

            % find uniform step
            if(nargin<2)
                maxDim = max(y1Max-y1Min, y2Max-y2Min);
                NPts = 10;
                dy = maxDim/NPts;
                dy = [dy;dy];
            elseif(length(dy)==1)
                dy = [dy;dy];
            end

            % find bounds for uniform grid
            y1Min0 = floor(y1Min/dy(1))*dy(1);
            y2Min0 = floor(y2Min/dy(2))*dy(2);
            y1Max0 = ceil(y1Max/dy(1))*dy(1);
            y2Max0 = ceil(y2Max/dy(2))*dy(2);

            % construct uniform grid
            y1 = (y1Min0:dy(1):y1Max0)';
            y2 = (y2Min0:dy(2):y2Max0)';

            N1u = length(y1);
            N2u = length(y2);
            y1_kv = kron(y1,ones(N2u,1));
            y2_kv = kron(ones(N1u,1),y2);

            %[y1_kv,y2_kv] = meshgrid(y1,y2);
            [x1_kv,x2_kv] = this.CompSpace(y1_kv,y2_kv);

            % remove those outside the shape (i.e. which get mapped outside
            % [-1,1] x [-1,1] in the computational space
            outOfBounds = (x1_kv > 1) | (x1_kv < -1) | (x2_kv > 1) | (x2_kv < -1);
            y1_kv(outOfBounds) = [];
            y2_kv(outOfBounds) = [];
            x1_kv(outOfBounds) = [];
            x2_kv(outOfBounds) = [];

            this.Pts.Uniform = struct('y1_kv',y1_kv,'y2_kv',y2_kv,'x1_kv',x1_kv,'x2_kv',x2_kv);

            UniformPts = this.Pts.Uniform;
            if(strcmp(this.polar,'polar'))
                this.Pts.y1_kv = tempPts.y1_kv;
                this.Pts.y2_kv = tempPts.y2_kv;
            end

        end


        %*************************************************************
        %******** Compute Interpolation Matrices *********************
        %*************************************************************

        function IP = InterpolationMatrix_Pointwise(this,y1P,y2P)
            IP = zeros(length(y1P),this.N1*this.N2);
            for i =1:length(y1P)
                [x1,x2] = CompSpace(this,y1P(i),y2P(i));
                h       = ComputeInterpolationMatrix(this,x1,x2,true,false);
                IP(i,:) = h.InterPol;
            end
        end

        % copy of function with different name for compatibility
        function IP = ComputeInterpolationMatrixPhys(this,y1P,y2P)
            IP = InterpolationMatrix_Pointwise(this,y1P,y2P);
        end

        function IP = SubShapePts(this,a_shapePts)
            IP = InterpolationMatrix_Pointwise(this,a_shapePts.y1_kv,...
                a_shapePts.y2_kv);
        end

        function IP = SubShapePtsCart(this,a_shapePts)
            pts = GetInvCartPts(this,a_shapePts.y1_kv,a_shapePts.y2_kv);
            IP = InterpolationMatrix_Pointwise(this,pts.y1_kv,pts.y2_kv);
        end

        function IP = InterpolationPlotCart(this,PlotArea,saveBool)

            y1C = PlotArea.y1Min + (PlotArea.y1Max-PlotArea.y1Min)*(0:(PlotArea.N1-1))'/(PlotArea.N1-1);
            y2C = PlotArea.y2Min + (PlotArea.y2Max-PlotArea.y2Min)*(0:(PlotArea.N2-1))'/(PlotArea.N2-1);

            y1C_kv = kronecker(y1C,ones(PlotArea.N2,1));
            y2C_kv = kronecker(ones(PlotArea.N1,1),y2C);

            pts         = GetInvCartPts(this,y1C_kv,y2C_kv);

            IP.InterPol = InterpolationMatrix_Pointwise(this,pts.y1_kv,pts.y2_kv);
            IP.pts1     = pts.y1_kv;
            IP.pts2     = pts.y2_kv;
            IP.Nplot1   = PlotArea.N1;
            IP.Nplot2   = PlotArea.N2;

            if((nargin == 3) && saveBool)
                this.Interp = IP;
            end

        end

        function IP = InterpolationPlot(this,PlotArea,saveBool)

            [x1min,x2min] = CompSpace(this,PlotArea.y1Min,PlotArea.y2Min);
            [x1max,x2max] = CompSpace(this,PlotArea.y1Max,PlotArea.y2Max);

            if(isfield(PlotArea,'N') && ~isfield(PlotArea,'N1'))
                PlotArea.N1 = PlotArea.N(1);
                PlotArea.N2 = PlotArea.N(2);
            end

            y1Plot = GetArray(x1min,x1max,PlotArea.N1);
            y2Plot = GetArray(x2min,x2max,PlotArea.N2);

            if(nargin == 3)
                IP = ComputeInterpolationMatrix(this,y1Plot,y2Plot,true,saveBool);
            else
                IP = ComputeInterpolationMatrix(this,y1Plot,y2Plot,true,false);
            end
        end

        function InterpolationPlotFlux(this,PlotArea)
            if(isfield(PlotArea,'NFlux'))
                N = PlotArea.NFlux;
            else
                N = 30;
                cprintf('m','!!PlotArea.NFlux is not set!!\n 30 is used as default in function "InterpolationPlotFlux" of "Shape" class.\n');
            end

            dy1     = (PlotArea.y1Max-PlotArea.y1Min);
            dy2     = (PlotArea.y2Max-PlotArea.y2Min);

            if(length(N)==1)
                deltaY1 = max(dy1,dy2)/N;
                deltaY2 = max(dy1,dy2)/N;
            else
                deltaY1 = dy1/N(1);
                deltaY2 = dy2/N(2);
            end

            PlotArea.N1 = round(dy1/deltaY1);
            PlotArea.N2 = round(dy2/deltaY2);

            this.InterpFlux = InterpolationPlotCart(this,PlotArea);
        end

        function InterPolUniform = InterpolationPlotUniform(this,dy)

            if(~isfield(this.Pts,'Uniform'))
                if(nargin<2)
                    this.ComputeUniformPts;
                else
                    this.ComputeUniformPts(dy);
                end
            end

            this.Interp.InterPolUniform = this.InterpolationMatrix_Pointwise(this.Pts.Uniform.y1_kv,this.Pts.Uniform.y2_kv);

            InterPolUniform = this.Interp.InterPolUniform;

        end

        %********************************************************
        %****************** Plotting functions*******************
        %********************************************************
        function [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = plotStreamlines(this,flux,startMask1,startMask2,opts)
            mask    = ((this.Pts.y1_kv <= max(this.Interp.pts1)) & ...
                (this.Pts.y1_kv >= min(this.Interp.pts1)) & ...
                (this.Pts.y2_kv <= max(this.Interp.pts2)) & ...
                (this.Pts.y2_kv >= min(this.Interp.pts2)));
            PtsYCart  = GetCartPts(this,this.Pts.y1_kv,this.Pts.y2_kv);

            fl_y1  = this.Interp.InterPol*flux(1:this.N1*this.N2);
            fl_y1  = reshape(fl_y1,this.Interp.Nplot2,this.Interp.Nplot1);

            fl_y2  = this.Interp.InterPol*flux(this.N1*this.N2+1:end);
            fl_y2  = reshape(fl_y2,this.Interp.Nplot2,this.Interp.Nplot1);

            yCart  = GetCartPts(this,this.Interp.pts1,this.Interp.pts2);
            y1M    = reshape(yCart.y1_kv,this.Interp.Nplot2,this.Interp.Nplot1);
            y2M    = reshape(yCart.y2_kv,this.Interp.Nplot2,this.Interp.Nplot1);

            if((nargin >= 4) && isfield(opts,'y2CartShift'))
                y2M = y2M + opts.y2CartShift;
            end

            if((nargin == 2) || isempty(startMask1))
                startMask1 = this.Ind.bound;
                startMask2 = this.Ind.bound;
            elseif((nargin == 3) || isempty(startMask2))
                startMask2 = startMask1;
            end

            if(islogical(startMask1))
                startMask1 = PtsYCart.y1_kv(mask&startMask1);
                startMask2 = PtsYCart.y2_kv(mask&startMask2);
            end

            h = streamline(y1M,y2M,fl_y1,fl_y2,startMask1,startMask2);
            if((nargin >= 5) && isfield(opts,'color'))
                set(h,'color',opts.color);
            end
            if((nargin >= 5) && isfield(opts,'linewidth'))
                set(h,'linewidth',opts.linewidth);
            end
        end

        function [y1_s,y2_s,fl_y1,fl_y2,hq] = plotFlux(this,flux,maskAdd,fl_norm,lw,c,plain,y2CartShift)

            flux = flux(:);

            nSpecies=size(flux,2);
            if(nSpecies == 1)
                nCol = 1;
            else
                nCol = 2;
            end

            nRows  = ceil(nSpecies/nCol);

            if(isempty(this.InterpFlux))
                pts = this.Pts;
                IP  = diag(length(this.Pts.y1_kv));
                mask    = ((this.Pts.y1_kv <= max(this.Interp.pts1)) & ...
                    (this.Pts.y1_kv >= min(this.Interp.pts1)) & ...
                    (this.Pts.y2_kv <= max(this.Interp.pts2)) & ...
                    (this.Pts.y2_kv >= min(this.Interp.pts2)));
                if(exist('maskAdd','var') && ~isempty(maskAdd))
                    mask   = (maskAdd & mask);
                end
            else
                pts.y1_kv = this.InterpFlux.pts1;
                pts.y2_kv = this.InterpFlux.pts2;
                IP        = this.InterpFlux.InterPol;
                mask      = true(size(this.InterpFlux.pts1));
            end

            if(sum(mask)==0)
                return;
            end

            yCart  = GetCartPts(this,pts.y1_kv,pts.y2_kv);
            y1_s   = yCart.y1_kv(mask);
            y2_s   = yCart.y2_kv(mask);

            if(nargin >= 8)
                y2_s = y2_s + y2CartShift;
            end

            xl    = [(min(y1_s)-0.5) (max(y1_s)+0.5)];
            yl    = [(min(y2_s)-0.5) (max(y2_s)+0.5)];

            %             if(exist('fl_norm','var') && ~isempty(fl_norm))
            %                 y1_s  = [NaN;y1_s];
            %                 y2_s  = [NaN;y2_s];
            %             end

            for iSpecies=1:nSpecies

                fl_y1     = IP*flux(1:this.N1*this.N2,iSpecies);
                fl_y2     = IP*flux(this.N1*this.N2+1:end,iSpecies);

                fl_y1 = fl_y1(mask,:);  fl_y2 = fl_y2(mask,:);

                %exclude close to zero fluxes
                abs_fl    = sqrt(fl_y1.^2 + fl_y2.^2);
                mark_zero = (abs_fl > max(0.01,0.2*max(abs_fl)));
                fl_y1 = fl_y1(mark_zero,:);
                fl_y2 = fl_y2(mark_zero,:);

                %                 if(exist('fl_norm','var') && ~isempty(fl_norm))
                %                     O = ones(1,size(fl_y1,2));
                %                     fl_y1 = [fl_norm*O;fl_y1];
                %                     fl_y2 = [0*O;fl_y2];
                %
                %                     y1_ss  = y1_s([true;mark_zero]);
                %                     y2_ss  = y2_s([true;mark_zero]);
                %                 else
                y1_ss  = y1_s([mark_zero]);
                y2_ss  = y2_s([mark_zero]);
                %                 end

                if(exist('fl_norm','var') && ~isempty(fl_norm))
                    fluxScale = fl_norm;
                else
                    fluxScale = 1;
                end

                if(exist('lw','var') && exist('c','var') && ~isempty(lw))
                    if(nSpecies == 1)
                        hq = quiver(y1_ss,y2_ss,fl_y1,fl_y2,fluxScale,'LineWidth',lw,'Color',c);
                    else
                        hq(iSpecies) = quiver(y1_ss,y2_ss,fl_y1,fl_y2,fluxScale,'LineWidth',lw,'Color',c{iSpecies});
                    end
                else
                    hq(iSpecies) = quiver(y1_ss,y2_ss,fl_y1,fl_y2,fluxScale);
                end

                hold on;

                if(nargin < 7 || ~strcmp(plain,'plain'))
                    xlim(xl); ylim(yl);
                    h = xlabel('$y_1$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
                    h = ylabel('$y_2$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);
                    pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);
                    if(nSpecies > 1)
                        title(['Species ' num2str(iSpecies)]);
                    end
                    set(gca,'fontsize',20);
                    set(gca,'linewidth',1.5);
                end
            end
            hold off;

        end

        function [y1M,y2M,VIM] = plot(this,V,options,optDetails)


            if((nargin >= 3) && ~(iscell(options)))
                options = {options};
            else
                options = {};
            end

            if(nargin<4)
                optDetails = struct;
            end

            %options: 'SC' , 'contour'
            %V = V(:);

            nSpecies=size(V,2);
            if(nSpecies == 1)
                nCol = 1;
            else
                nCol = 2;
            end

            if(IsOption(options,'flux'))
                Interp = this.InterpFlux;
            else
                Interp = this.Interp;
            end

            nRows  = ceil(nSpecies/nCol);
            yCart  = GetCartPts(this,Interp.pts1,Interp.pts2);


            if((nargin >= 4) && isfield(optDetails,'y2CartShift'))
                yCart.y2_kv = yCart.y2_kv + optDetails.y2CartShift;
            end

            yCart1 = yCart.y1_kv;
            yCart2 = yCart.y2_kv;
            yCart1 = yCart1(isfinite(yCart1));
            yCart2 = yCart2(isfinite(yCart2));

            xl     = [min(yCart1) max(yCart1)];
            yl     = [min(yCart2) max(yCart2)];

            if((nargin >= 3) && IsOption(options,'comp'))
                [x1,x2]  = CompSpace(this,Interp.pts1,Interp.pts2);
                y1M    = reshape(x1,Interp.Nplot2,Interp.Nplot1);
                y2M    = reshape(x2,Interp.Nplot2,Interp.Nplot1);
            else
                y1M    = reshape(yCart.y1_kv,Interp.Nplot2,Interp.Nplot1);
                y2M    = reshape(yCart.y2_kv,Interp.Nplot2,Interp.Nplot1);
            end

            if((nargin >= 4) && isfield(optDetails,'linewidth'))
                lw = optDetails.linewidth;
            else
                lw = 1.0;
            end
            for iSpecies=1:nSpecies

                if(nSpecies > 1)
                    subplot(nRows,nCol,iSpecies);
                end

                if((nargin >= 3) && IsOption(options,'flux')) %if its a flux
                    VI = blkdiag(Interp.InterPol,Interp.InterPol)*V(:,iSpecies);
                elseif( (size(V,1) == length(Interp.pts1)) && (length(Interp.pts1) ~= length(this.Pts.y1_kv)) )
                    VI = V(:,iSpecies);
                else
                    VI = Interp.InterPol*V(:,iSpecies);
                end

                if((nargin >= 3) && IsOption(options,'contour'))
                    if(nargin >= 4)
                        VIM   = reshape(VI,Interp.Nplot2,Interp.Nplot1);
                        if(isfield(optDetails,'nContours'))
                            ns = optDetails.nContours;
                            if(length(optDetails.nContours)==1)
                                ns = ns*[1,1];
                            end
                            [C,h] = contour(y1M,y2M,VIM,ns,'linewidth',lw);
                        else
                            [C,h] = contour(y1M,y2M,VIM,...
                                'linewidth',lw);
                        end
                        if(isfield(optDetails,'linecolor'))
                            set(h,'Color',optDetails.linecolor);
                        end
                        if(isfield(optDetails,'clabel') && optDetails.clabel)
                            hc = clabel(C,h); set(hc,'fontsize',15);
                        end

                    else
                        VIM   = reshape(VI,Interp.Nplot2,Interp.Nplot1);
                        [C,h] = contour(y1M,y2M,VIM,8,...
                            'linewidth',lw);
                    end
                    % hh = clabel(C,h,'fontsize',20);
                elseif((nargin >= 3) && IsOption(options,'color'))
                    VIM   = reshape(VI,Interp.Nplot2,Interp.Nplot1);
                    pcolor(y1M,y2M,VIM);
                    shading flat;
                    %                     colorbar;
                    %                     minV = min(min(VIM));
                    %                     maxV = max(max(VIM));
                    %                     if(abs(maxV - minV)> 0)
                    %                        set(gca, 'CLim', [minV maxV]);
                    %                    end
                elseif((nargin >= 3) && IsOption(options,'flux'))
                    VIM.fl_y1     = VI(1:end/2,iSpecies);
                    VIM.fl_y2     = VI(end/2+1:end,iSpecies);

                    if((nargin >= 4) && isfield(optDetails,'linecolor'))
                        c = optDetails.linecolor;
                    else
                        c = 'k';
                    end
                    y1M = yCart.y1_kv;
                    y2M = yCart.y2_kv;
                    quiver(y1M,y2M,VIM.fl_y1,VIM.fl_y2,'LineWidth',lw,'Color',c);
                else
                    VIM   = reshape(VI,Interp.Nplot2,Interp.Nplot1);
                    h = surf(y1M,y2M,VIM); %mesh

                    if(nargin>3)
                        if(isfield(optDetails,'linecolor'))
                            set(h,'FaceColor',optDetails.linecolor);
                            set(h,'FaceAlpha',0.5);
                        end
                        if(isfield(optDetails,'edgecolor'))
                            set(h,'EdgeColor',optDetails.edgecolor);
                        end
                        if(isfield(optDetails,'zLim'))
                            set(gca,'zLim',optDetails.zLim);
                            set(gca,'cLim',optDetails.zLim);
                        end
                    end
                end
                hold on;

                if((nargin >= 3) && IsOption(options,'SC'))
                    hold on;

                    if(IsOption(options,'comp'))
                        pts.y1_kv = this.Pts.x1_kv;
                        pts.y2_kv = this.Pts.x2_kv;
                        mask      = true(size(pts.y1_kv));
                    else
                        pts = GetCartPts(this);
                        mask    = ((pts.y1_kv <= xl(2)) & ...
                            (pts.y1_kv >= xl(1)) & ...
                            (pts.y2_kv <= yl(2)) & ...
                            (pts.y2_kv >= yl(1)));
                    end
                    Vp      = V(:,iSpecies);
                    h       = scatter3(pts.y1_kv(mask),pts.y2_kv(mask),...
                        Vp(mask),'o');
                    set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
                end

                if((nargin < 4) || ~isfield(optDetails,'reshape') || optDetails.reshape)
                    if((nargin >= 3) && IsOption(options,'comp'))
                        xlim([-1,1]); ylim([-1,1]);
                    else
                        xlim(xl); ylim(yl);
                    end
                    if(nargin>=4 && isfield(optDetails,'sigma') && optDetails.sigma)
                        xlabel('$x/\sigma$','Interpreter','Latex','fontsize',20);
                        ylabel('$y/\sigma$','Interpreter','Latex','fontsize',20);
                    elseif(nargin>=4 && isfield(optDetails,'xLabel'))
                        xlabel(optDetails.xLabel,'Interpreter','Latex','fontsize',20);
                        ylabel(optDetails.yLabel,'Interpreter','Latex','fontsize',20);
                    else
                        xlabel('$x$','Interpreter','Latex','fontsize',20);
                        ylabel('$y$','Interpreter','Latex','fontsize',20);
                    end
                    pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);
                    if(nargin>=4 && isfield(optDetails,'title'))
                        title(optDetails.title,'Interpreter','Latex','fontsize',20);
                    else
                        if(nSpecies > 1)
                            title(['Species ' num2str(iSpecies)]);
                        end
                    end

                    if(isfield(optDetails,'fontsize'))
                        set(gca,'fontsize',optDetails.fontsize);
                    else
                        set(gca,'fontsize',20);
                    end

                    set(gca,'linewidth',1.5);
                    hold off;
                end
            end

            if(nargout == 0)
                y1M = [];
            end

        end

        function h = PlotGrid(this,InitFig)

            ptsCart = GetCartPts(this);
            h = scatter(ptsCart.y1_kv,ptsCart.y2_kv,'ob','HandleVisibility','off');
            set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');

            if((nargin == 2) && InitFig)
                y1lim = [min(ptsCart.y1_kv) max(ptsCart.y1_kv)];
                y2lim = [min(ptsCart.y2_kv) max(ptsCart.y2_kv)];

                xlim(y1lim + [-0.5 0.5]);
                ylim(y2lim + [-0.5 0.5]);

                pbaspect([(y1lim(2)-y1lim(1)) (y2lim(2)-y2lim(1)) 1]);
                set(gca,'fontsize',15);
            end
        end

        function PlotNormals(this)
            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;

            dir1 = [ones(size(y1)); zeros(size(y2))];
            dir2 = [zeros(size(y1)); ones(size(y2))];

            if(isfield(this.Ind, 'normalCart'))
                normal = this.Ind.normalCart;
            else
                normal = this.Ind.normal;
            end

            bound = this.Ind.bound;

            if(strcmp(this.polar,'polar'))
                ptsCart = this.GetCartPts(y1,y2);
                y1 = ptsCart.y1_kv;
                y2 = ptsCart.y2_kv;
            end

            v1 = normal * dir1;
            v2 = normal * dir2;

            vnorm = sqrt(v1.^2 + v2.^2);
            v1 = v1./vnorm;
            v2 = v2./vnorm;

            y1 = y1(bound);
            y2 = y2(bound);

            quiver(y1, y2, v1, v2, 'AutoScale','off');

            hold on
        end

        function PlotBound(this,boundList)
            if(nargin==1)
                boundList = {'bound'};
            end

            colourList = {'r','b','m','c'};

            %             y1 = this.Pts.y1_kv;
            %             y2 = this.Pts.y2_kv;

            y1 = this.Pts.y1_kv;
            y2 = this.Pts.y2_kv;

            if(strcmp(this.polar,'polar'))
                ptsCart = this.GetCartPts(y1,y2);
                y1 = ptsCart.y1_kv;
                y2 = ptsCart.y2_kv;
            end

            nFields = length(boundList);

            for iField = 1:nFields
                mask = this.Ind.(boundList{iField});
                p1 = scatter(y1(mask),y2(mask),'o','MarkerEdgeColor',colourList{iField},'MarkerFaceColor',colourList{iField});
                hold on
            end
            legend(p1,'bound');
        end

        function PlotIsoline(this,x,y1y2)
            %xI = (-1:0.01:1)';
            if(strcmp(y1y2,'y1'))
                xI = (min(this.Pts.x1):0.01:1)';
                [y1_kv,y2_kv] = PhysSpace(this,x*ones(size(xI)),xI);
            elseif(strcmp(y1y2,'y2'))
                xI = (min(this.Pts.x2):0.01:1)';
                [y1_kv,y2_kv] = PhysSpace(this,xI,x*ones(size(xI)));
            else
                return;
            end
            ptsC = GetCartPts(this,y1_kv,y2_kv);
            plot(ptsC.y1_kv,ptsC.y2_kv,'r','linewidth',1.5); hold on;
        end

        function PlotGridLines(this,optDetails)

            if((nargin >= 2) && isfield(optDetails,'y2CartShift'))
                y2CartShift = optDetails.y2CartShift;
            else
                y2CartShift = 0;
            end

            if((nargin >= 2) && isfield(optDetails,'nthGridLines'))
                nthGridLines = optDetails.nthGridLines;
            else
                nthGridLines = 1;
            end

            xI = (-1:0.05:1)';
            x1I = (min(this.Pts.x1):0.05:1)';
            x2I = (min(this.Pts.x2):0.05:1)';

            O  = ones(size(xI));
            O1  = ones(size(x1I));
            O2  = ones(size(x2I));
            %(1) Plot x1-isolines
            for i1=1:nthGridLines:this.N1
                [y1_kv,y2_kv] = PhysSpace(this,this.Pts.x1(i1)*O2,x2I);
                GL_CartPts = GetCartPts(this,y1_kv,y2_kv);
                plot(GL_CartPts.y1_kv,GL_CartPts.y2_kv+y2CartShift,'b'); hold on;
            end

            %(2) Plot x2-isolines
            for i2=1:nthGridLines:this.N2
                [y1_kv,y2_kv] = PhysSpace(this,x1I,this.Pts.x2(i2)*O1);
                GL_CartPts    = GetCartPts(this,y1_kv,y2_kv);
                plot(GL_CartPts.y1_kv,GL_CartPts.y2_kv+y2CartShift,'b'); hold on;                %'linewidth',1.
            end

        end

        function [I,w,weights,IP,pts] = doIntFLine(this,y1P,y2P,f,TRAP_CHEB)

            N         = 500;
            if((nargin < 5) || (strcmp(TRAP_CHEB,'TRAP')))
                xi        = (0:N)'/N;
                weightX   = [0.5,ones(1,length(xi)-2),0.5]/N;
            elseif(strcmp(TRAP_CHEB,'CHEB'))
                [xi,weightX] = ClenCurtFlip(N);
                xi      = (1+xi)/2;
                weightX = weightX/2;
            end

            if( (y1P(1) == inf) && (y1P(2) == inf))
                d1 = 0;
            else
                d1 = y1P(2)-y1P(1);
            end

            if( (y2P(1) == inf) && (y2P(2) == inf))
                d2 = 0;
            else
                d2 = y2P(2)-y2P(1);
            end


            pts.y1_kv = y1P(1) + xi*d1;
            pts.y2_kv = y2P(1) + xi*d2;

            IP        = SubShapePtsCart(this,pts);

            dist      = sqrt(d1^2 + d2^2);
            weights   = weightX*dist;

            w         = weights*IP;

            if(isempty(f))
                I = 0;
            else
                I = w*f;
            end
        end

        function [f_p,pts] = plotLine(this,y1P,y2P,f,opts)

            plain = false;
            color = 'k';
            dist0 = false;
            CART  = true;

            if((nargin > 4))
                if(isfield(opts,'plain'))
                    plain = opts.plain;
                end

                if(isfield(opts,'color'))
                    color = opts.color;
                end

                if(isfield(opts,'dist0'))
                    dist0 = opts.dist0;
                end

                if(isfield(opts,'CART'))
                    CART = opts.CART;
                end
            end

            xi        = (0:0.002:1)';

            if(y1P(2) == y1P(1))
                dy1 = 0;
            else
                dy1 = y1P(2) - y1P(1);
            end
            pts.y1_kv = y1P(1) + xi*dy1;

            if(y2P(2) == y2P(1))
                dy2 = 0;
            else
                dy2 = y2P(2) - y2P(1);
            end

            pts.y2_kv = y2P(1) + xi*dy2;

            if((y1P(1) == y1P(2)))
                dist = sqrt( (xi*(y2P(2)-y2P(1))).^2);
            elseif((y2P(1) == y2P(2)))
                dist = sqrt( (xi*(y1P(2)-y1P(1))).^2);
            else
                dist = sqrt( (xi*(y1P(2)-y1P(1))).^2 + ...
                    (xi*(y2P(2)-y2P(1))).^2 );
            end
            pts.z = dist;

            %Plot Grid Lines that are crossed:
            %(a) Get x-Values of initial points
            if(CART)
                ptY_start = GetInvCartPts(this,y1P(1),y2P(1));
                ptY_end   = GetInvCartPts(this,y1P(2),y2P(2));
            else
                ptY_start.y1_kv = y1P(1);
                ptY_start.y2_kv = y2P(1);

                ptY_end.y1_kv = y1P(2);
                ptY_end.y2_kv = y2P(2);
            end

            %(b) Find y1-points which lie between start and end point
            markG1 = ( (this.Pts.y1 > ptY_start.y1_kv) & (this.Pts.y1 < ptY_end.y1_kv));
            xmG1   = (this.Pts.y1(markG1) - ptY_start.y1_kv)/(ptY_end.y1_kv-ptY_start.y1_kv);
            ptsG1.y1_kv   = y1P(1) + xmG1*(y1P(2)-y1P(1));
            ptsG1.y2_kv   = y2P(1) + xmG1*(y2P(2)-y2P(1));

            distG1        = sqrt( (xmG1*(y1P(2)-y1P(1))).^2 + ...
                (xmG1*(y2P(2)-y2P(1))).^2 );

            %(c) Find y2-points which lie between start and end point
            markG2 = ( (this.Pts.y2 > ptY_start.y2_kv) & (this.Pts.y2 < ptY_end.y2_kv));
            xmG2   = (this.Pts.y2(markG2) - ptY_start.y2_kv)/(ptY_end.y2_kv-ptY_start.y2_kv);
            ptsG2.y1_kv   = y1P(1) + xmG2*(y1P(2)-y1P(1));
            ptsG2.y2_kv   = y2P(1) + xmG2*(y2P(2)-y2P(1));

            distG2        = sqrt( (xmG2*(y1P(2)-y1P(1))).^2 + ...
                (xmG2*(y2P(2)-y2P(1))).^2 );

            if(CART)
                IP       = SubShapePtsCart(this,pts);
                IPG1     = SubShapePtsCart(this,ptsG1);
                IPG2     = SubShapePtsCart(this,ptsG2);
            else
                IP       = SubShapePts(this,pts);
                IPG1     = SubShapePts(this,ptsG1);
                IPG2     = SubShapePts(this,ptsG2);
            end

            if((y1P(1) == y1P(2)) && ~dist0)
                offset = y2P(1);
                xLab_Txt = 'y_2';
            elseif(y2P(1) == y2P(2)  && ~dist0)
                offset   = y1P(1);
                xLab_Txt = 'y_1';
            else
                offset   = 0;
                xLab_Txt = 'Distance from point y0';
            end

            f_p = IP*f;
            plot(offset + dist,f_p,'color',color); hold on;
            xlabel(xLab_Txt);
            if(~plain)
                plot(offset + distG1,IPG1*f,'o','MarkerEdgeColor','k','MarkerFaceColor','g');
                plot(offset + distG2,IPG2*f,'o','MarkerEdgeColor','k','MarkerFaceColor','g');
                ptsStr = ['(',num2str(y1P(1)),',',num2str(y2P(1)),') to (',num2str(y1P(2)),',',num2str(y2P(2)),')'];
                title(['Values on line from ',ptsStr]);
            end
        end

        function makeMovie(this,rhot,options,optDetails)
            if(nargin<3)
                options = {};
            end
            if(nargin<4)
                optDetails = struct;
            end

            if(isfield(optDetails,'doPdfs'))
                doPdfs = optDetails.doPdfs;
            else
                doPdfs = false;
            end

            if(doPdfs)
                if(isfield(optDetails,'pdfDir'))
                    pdfDir = optDetails.pdfDir;
                else
                    pdfDir = ['tempMovie' filesep];
                end
                [nd,pdfFileNames] = setupPdfMovie(rhot,pdfDir);
            end

            if(isfield(optDetails,'doMp4'))
                doMp4 = optDetails.doMp4;
            else
                doMp4 = false;
            end

            % open an mp4 file to write to
            if(doMp4)
                if(isfield(optDetails,'mp4File'))
                    v = VideoWriter(optDetails.mp4File,'MPEG-4');
                else
                    v = VideoWriter('movieDemo.mp4','MPEG-4');
                end
                v.FrameRate=10;
                v.Quality = 100;
                open(v)
            end


            if(isfield(optDetails,'addTimes'))
                addTimes = optDetails.addTimes;
                plotTimes = optDetails.plotTimes;
                if(isfield(optDetails,'title'))
                    baseTitle = optDetails.title;
                else
                    baseTitle = '';
                end
            else
                addTimes = false;
            end


            nTimes = size(rhot,2);
            for iT = 1:nTimes

                if(addTimes)
                    tString = sprintf('%0.1f',plotTimes(iT));
                    if(isempty(baseTitle))
                        optDetails.title = ['$t =' tString '$'];
                    else
                        optDetails.title = [baseTitle ', $t =' tString '$'];
                    end
                end

                this.plot(rhot(:,iT),options,optDetails);

                % save the mp4 frame
                if(doMp4)
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end

                % save the pdfs
                if(doPdfs)
                    pdfFileNames = savePdfFrame(pdfDir,pdfFileNames,iT,nd,gcf);
                else
                    pause(0.1);
                end

            end

            % close & write the mp4 file
            if(doMp4)
                close(v);
            end

            % make the full movie pdf
            if(doPdfs)
                makePdfMovie(pdfDir,pdfFileNames);
            end

        end

        function [hf,ha,plotFile] = DocumentationPlot(this,opts)
            hf = []; % by default we don't have a plot
            ha = []; % by default we don't have a plot
            plotFile = []; % by default we don't have a plot
        end

    end    

end
