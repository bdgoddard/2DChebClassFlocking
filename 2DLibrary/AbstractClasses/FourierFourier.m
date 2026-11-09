classdef (Abstract) FourierFourier < Shape
    
    properties (Access = public)
        FFTMatrix
        IFFTMatrix
    end
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = FourierFourier(N1,N2)
            this@Shape(N1,N2);                                    
             
            if(mod(N1,2) == 1)
                exception = MException('FourierFourierClass:Constructor','N1 has to be even');
                throw(exception);                 
            end                          
            if(mod(N2,2) == 1)
                exception = MException('FourierFourierClass:Constructor','N2 has to be even');
                throw(exception);                 
            end                          

            this.Pts.x1 = (0:(N1-1))'/(N1);   
            this.Pts.x2 = (0:(N2-1))'/(N2);   

            this.Pts.N1 = N1;
            this.Pts.N2 = N2;

            [this.FFTMatrix,this.IFFTMatrix] = getFFTMatrices(this,N1,N2);     
        end      
    end
    
    %**********************************************
    %************   Initializations    *************
    %**********************************************
    methods (Access = public) 
        function [int,area] = ComputeIntegrationVector(this)
            [h1,dy1]        = this.PhysSpace1(this.Pts.x1);
            [h1,dy2]        = this.PhysSpace2(this.Pts.x2);
            %1/Nj weight of integration in jth direction
            this.Int = kron((dy1').*ones(1,this.N1)/this.N1,(dy2').*ones(1,this.N2)/this.N2);            
            int = this.Int;
        end        

        function Diff = ComputeDifferentiationMatrix(this,pts)
            
            if(nargin == 1)
                pts = this.Pts;                
            end
            
            N1 = length(pts.x1);
            N2 = length(pts.x2);

            %Differentiation matrix in x1
            h         = 2*pi/N1;       
            column    = [0 .5*(-1).^(1:N1-1).*cot((1:N1-1)*h/2)];
            Diff1.Dx  = toeplitz(column,column([1 N1:-1:2]))*2*pi;% 2pi is needed to go back to the -1;1 - domain        

            h         = 2*pi/N1;       
            column2   = [(-(pi^2)/(3*h^2) - 1/6)  -.5*(-1).^(1:N1-1)./((sin((1:N1-1)*h/2)).^2)];
            Diff1.DDx = toeplitz(column2,column2([1 N1:-1:2]))*(2*pi)^2;        
            
            %Differentiation matrix in x2
            h         = 2*pi/N2;       
            column    = [0 .5*(-1).^(1:N2-1).*cot((1:N2-1)*h/2)];
            Diff2.Dx  = toeplitz(column,column([1 N2:-1:2]))*2*pi;% 2pi is needed to go back to the -1;1 - domain        

            h         = 2*pi/N2;       
            column2   = [(-(pi^2)/(3*h^2) - 1/6)  -.5*(-1).^(1:N2-1)./((sin((1:N2-1)*h/2)).^2)];
            Diff2.DDx = toeplitz(column2,column2([1 N2:-1:2]))*(2*pi)^2;        


            Sel = {'Dy1' ;'DDy1' ; 'Dy2'; 'DDy2' ;...
                  'Dy1Dy2'; 'Lap' ;'grad' ;'div';...
                  'gradDiv'; 'LapVec'}; 

            this.Diff = PhysicalDerivatives(this,pts,Sel,Diff1,Diff2,2);
            Diff = this.Diff;
            %Dx1,Dx2,DDx1,DDx2);        

        end  

        function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)

            N1 = this.N1;
            N2 = this.N2;
            
            %1st variable: Fourier
            LL    = 2*pi;
            M1     = floor((N1+1)/2);

            %%FFT Interpolation
            n11   = 1:(M1-1);       % first half of points, excluding zero
            n12   = (M1+1):(N1-1);  % second half

            % This uses that NT+1 is odd (ie NT is even)
            phi      = interp1*2*pi;            
            Interp1  = [exp(2*pi*1i*phi/LL*0) , exp(2*pi*1i*phi/LL*n11) , cos(2*pi*M1*phi/LL) , exp(2*pi*1i*phi/LL*(n12-N1))]/N1;

            
            %2nd variable: Fourier
            M2     = floor((N2+1)/2);

            %%FFT Interpolation
            n21   = 1:(M2-1);       % first half of points, excluding zero
            n22   = (M2+1):(N2-1);  % second half

            % This uses that NT+1 is odd (ie NT is even)
            phi      = interp2*2*pi;            
            Interp2  = [exp(2*pi*1i*phi/LL*0) , exp(2*pi*1i*phi/LL*n21) , cos(2*pi*M2*phi/LL) , exp(2*pi*1i*phi/LL*(n22-N2))]/N2;

            InterPol = real(kron(Interp1,Interp2)*this.FFTMatrix);
            
            interp_y1 = this.PhysSpace1(interp1);
            interp_y2 = this.PhysSpace2(interp2);
            % kron form of the two interpolations                                                            
            Interp = struct('InterPol',InterPol,...
                            'interp1',interp_y1,'interp2',interp_y2,...             
                            'Nplot1',length(interp1),...
                            'Nplot2',length(interp2),...
                            'N1',this.N1,'N2',this.N2);                   
                        
            if((nargin >= 4) && fullInterpGrid)
                Interp.pts1 = kron(interp_y1,ones(size(interp_y2)));
                Interp.pts2 = kron(ones(size(interp_y1)),interp_y2);             
            end
            if((nargin >= 5) && saveBool)
                this.Interp = Interp;
            end                        
        end        
        
        function Ind = ComputeIndices(this)
            Ind      = GetIndicesFourierFourier(this);
            this.Ind = Ind;
        end    

        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)

            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            N1  = this.N1;  
            N2  = this.N2;
            Pts = this.Pts;
            Int = this.Int;  % 1 x N1*N2

            if(useDistance)
                fPTemp = f(GetDistance(this,Pts.y1_kv,Pts.y2_kv));
            else
                fPTemp = f(Pts.y1_kv,Pts.y2_kv);
            end
            
            fDim = size(fPTemp);
            
            nElts = prod(fDim(2:end));
            
            IntT = Int.';  % N1*N2 x 1
            
            IntT = IntT(:,ones(1,nElts)); % N1*N2 x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([N1*N2,N1*N2,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:(N1*N2) 
                if(useDistance)
                    fP          = f(GetDistance(this,Pts.y1_kv(i) - Pts.y1_kv,Pts.y2_kv(i) - Pts.y2_kv));
                else
                    fP          = f(GetSignedDistance1(this,Pts.y1_kv(i) - Pts.y1_kv),GetSignedDistance2(this,Pts.y2_kv(i) - Pts.y2_kv));
                end
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;
            
            if((nargin >= 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end

        end        
        
        function [M_ConvConv,C1,C2,C1Line,C2Line] = ComputeConvolutionConvolutionMatrix(this,f1,f2,shapeParams1,shapeParams2)
            
            % Convolution in y1 (Fourier) potentially with cutoff
            
            % set up Fourier line
            geom1.N = this.N1;
            geom1.yMin = this.y1Min;
            geom1.yMax = this.y1Max;
            f1Line = FourierLine(geom1);

            % compute convolution in y1
            f1Line.ComputeIntegrationVector;
            C1Line = f1Line.ComputeConvolutionMatrix(f1,shapeParams1);
            C1 = sparse(kron(C1Line,eye(this.N2)));
            
            % Convolution in y2 (Fourier) potentially with cutoff
            
            % set up Fourier line
            geom2.N = this.N2;
            geom2.yMin = this.y2Min;
            geom2.yMax = this.y2Max;
            f2Line = FourierLine(geom2);

            % compute convolution (which requires the integration vector)
            f2Line.ComputeIntegrationVector;
            C2Line = f2Line.ComputeConvolutionMatrix(f2,shapeParams2);
            C2 = sparse(kron(eye(this.N1),C2Line));
            
            % Full convolution
            M_ConvConv = C1*C2;
        end

        function M_conv = ComputeConvolutionMatrix_Pointwise(this,f,shapeParams)
            % Convolution \int_D(x) f(x-y)g(y) dy
            % Here D(x) is the intersection of the subshape defined in
            % shapeParams with the periodic box, i.e., just the whole shape
            
            disp('Computing Convolution matrices...');
    
            % determine if the kernel uses a distance or the two
            % coordinates
            if(nargin(f)==1)
                useDistance = true;
            else
                useDistance = false;
            end

            if(isfield(shapeParams,'shape'))
                % determine the shape we're using
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
            
            % allow computing for vector functions
            fDim = size(fPTemp,2);
            
            N1  = this.N1;  N2  = this.N2;
            N1N2 = N1*N2;
            M_conv = zeros(N1*N2,N1*N2,fDim);
            
            wb = waitbar(0,'Computing Convolution Matrices');
            
            % get gridpoints to centre the domain on
            y1 = Pts.y1_kv;
            y2 = Pts.y2_kv;
             
            % set up the subshape (note that we don't need to do the
            % intersection
            subShape = subshapeFn(shapeParams);
            % compute the integration weights
            subShape.ComputeIntegrationVector;

            % determine argument for f(x-y)
            % note that we will shift the points in subShape by y1(i) and
            % y2(i), which correspond to x1 and x2
            % we then compute x-y which is -subShapePts
            % so f doesn't depeond on the point x
            % this should be clear from the argument only depending on the
            % distance from the centre of the shape
            dy1 = - subShape.Pts.y1_kv;
            dy2 = - subShape.Pts.y2_kv;

            if(useDistance)
                fP          = f(GetDistance(this,dy1,dy2));
            else
                fP          = f(dy1,dy2);
            end


            for i=1:N1N2
                
                waitbar(i/N1N2,wb);

                % shift the origin to the current point to
                % determine the points in the intersected domain
                subConvShapePts = shiftPointsCart(subShape.Pts,y1(i),y2(i));

                % the interpolation onto D(x) from the full domain
                IP_i     = SubShapePts(this,subConvShapePts);
                
                % the integration weights for the subshape
                Int_i    = subShape.Int;
                
                % form the convolution matrix
                for iF = 1:fDim
                    M_conv(i,:,iF) = (Int_i.*fP(:,iF)')*IP_i;
                end

            end

            M_conv(isnan(M_conv)) = 0;
                
            disp('...done.');
            
            close(wb);
            
        end % convolution               



        % overridden from Shape
        function d  = GetDistance(this,pts_y1,pts_y2)             
            % pts_y1 and Pts_y are the y1 and y2 components of the
            % non-periodic distances
            
            if(nargin == 2)
                tempPts = pts_y1;
                pts_y1 = tempPts.y1_kv;
                pts_y2 = tempPts.y2_kv;
            end
            
            d1 = this.GetSignedDistance1(pts_y1);
            d2 = this.GetSignedDistance2(pts_y2);

            d       = sqrt(d1.^2 + d2.^2);
        end

        function d  = GetSignedDistance1(this,y1)             
            d = y1;
            mask = ( abs(d) > this.L1/2 );
            d(mask) = d(mask) - sign(d(mask))*this.L1;
        end
        function d  = GetSignedDistance2(this,y2)             
            d = y2;
            mask = ( abs(d) > this.L2/2 );
            d(mask) = d(mask) - sign(d(mask))*this.L2;            
        end
        
        function [FFTMatrix,IFFTMatrix] = getFFTMatrices(this,N1,N2)

            n1 = (1:N1) -1;
            k1 = (1:N1) -1;

            FFTMatrix1 = exp( -1i * 2*pi *n1'*k1 /N1);
            IFFTMatrix1 = 1/N1 * exp( 1i * 2*pi *n1'*k1 /N1);

            n2 = (1:N2) -1;
            k2 = (1:N2) -1;

            FFTMatrix2 = exp( -1i * 2*pi *n2'*k2 /N2);
            IFFTMatrix2 = 1/N2 * exp( 1i * 2*pi *n2'*k2 /N2);

            FFTMatrix = kron(FFTMatrix1,FFTMatrix2);
            IFFTMatrix = kron(IFFTMatrix1,IFFTMatrix2);

        end        
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
           
            [y1_kv,dy1] = PhysSpace1(this,x1);
            [y2_kv,dy2] = PhysSpace2(this,x2);
            
            if(nargout > 2)
                n        = length(x1);
                J        = zeros(n,2,2);
                J(:,1,1) = dy1;
                J(:,2,2) = dy2;                
            end
            
            if(nargout > 3)
                dH1 = 0;
                dH2 = 0;
                disp('SpectralFourierClass:Comp_to_Phys: dH1/2 NOT YET IMPLEMENTED');
            end
            
        end
        function [x1,x2] = CompSpace(this,y1,y2)
            x1 = CompSpace1(this,y1);
            x2 = CompSpace2(this,y2);
        end
    end
    
    %**********************************************
    %************   Mapping functions *************
    %**********************************************
    methods (Abstract = true,Access = public)         
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(x);
         [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(x);
         x = CompSpace1(y);
         x = CompSpace2(y);
    end  
    
end


