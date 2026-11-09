classdef (Abstract) Fourier < Interval
    %FOURIER Class to encapsulate Fourier grid. It derives from Interval
    %class so it will show all the expected functionality given by that
    %another abstract class.
    
    properties (Access = public)
        FFTMatrix
        IFFTMatrix
    end
    
    methods 
        function this = Fourier(N)
            
            if(mod(N,2) == 1)              
                msgID = 'FourierClass:Constructor';
                msgtext = 'Number of points, N, must be even.';
                ME = MException(msgID,msgtext,N); 
                throwAsCaller(ME);
            end
            
            this@Interval(N);
            this.Pts.x = (0:N-1)'/N; % computational domain
            %this.Pts.y = (2*pi)*this.Pts.x; % physical domain
            
            [this.FFTMatrix,this.IFFTMatrix] = this.getFFTMatrices;
        end      
    end
    
    methods (Access = public)
        
        function Diff = ComputeDifferentiationMatrix(this)        
            CompDiff = FourierDiff(this.Pts.x,2);             
            Diff = PhysicalDerivatives_1D(@this.PhysSpace,this.Pts.x,CompDiff);   
            
            this.Diff = Diff;              
        end 
        
        function Int =  ComputeIntegrationVector(this)
            [~,wInt] = FourierInt(this.N);
            [~,J] = PhysSpace(this,this.Pts.x);
            
            J(J==inf)  = 0;  
            J(J==-inf)  = 0;  
            J(isnan(J))= 0;  
            
            Int = wInt.*(J.');
            this.Int = Int;
        end
        
        function Interp = ComputeInterpolationMatrix(this,interp,saveBool)                       
            
            InterPol = fftInterpMatrix(this.Pts.x,interp);               
            Interp = struct('InterPol',InterPol,...
                            'Nplot',length(interp),...
                            'N',this.N,...
                            'pts',this.PhysSpace(interp),...
                            'xPlot',interp,...
                            'yPlot',this.PhysSpace(interp));

            if((nargin >= 3) && saveBool)
                this.Interp = Interp;
            end  
            
        end  
                
        function Ind    = ComputeIndices(this)
            Ind = GetIndicesFourier(this);
            this.Ind = Ind;
        end
        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            N  = this.N; 
            y = this.Pts.y;
            Int = this.Int;  % 1 x N

            % find size of function matrix by making a dummy version
            fPTemp = f(y);
            fDim   = size(fPTemp);
            nElts  = prod(fDim(2:end)); % first dimension stores values
            
            IntT = Int.';  % N x 1
            
            IntT = IntT(:,ones(1,nElts)); % N x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([N,N,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:N 
                fP          = f(GetDistance(this,y(i)-y));  
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            %M_conv(isnan(M_conv)) = 0;
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end                

        
        function M_conv = ComputeConvolutionMatrix_old(this,f,saveBool)
            %M_conv = [];
            
            N  = this.N; 
            y = this.Pts.y;
            Int = this.Int;
            
            M_conv = zeros(N);

            for i=1:N 
                M_conv(i,:)  = Int.*f(GetDistance(this,y(i)-y)).';
            end
            
            this.Conv = M_conv;
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end
        
        function d  = GetDistance(this,y)
            d = y;
            mask = ( abs(d) > this.L/2 );
            d(mask) = d(mask) - sign(d(mask))*this.L;
            %d = mod(y,this.L);
        end

        
        function [FFTMatrix,IFFTMatrix] = getFFTMatrices(this,saveBool)
            N = this.N;
            n = (0:N-1);
            k = (0:N-1);
            
            FFTMatrix = exp( -1i * 2*pi *n'*k /N);            
            IFFTMatrix = 1/N * exp( 1i * 2*pi *n'*k /N);
            
            if((nargin >= 2) && saveBool)
                this.FFTMatrix  = FFTMatrix;
                this.IFFTMatrix = IFFTMatrix;
            end  
        end
    end
end

