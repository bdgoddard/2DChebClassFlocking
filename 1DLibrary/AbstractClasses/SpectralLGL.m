classdef (Abstract) SpectralLGL < Interval
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = SpectralLGL(N)
             this@Interval(N);             
             [x,~,P] = LGLPointsWeights(N-1); 
             this.Pts.x = x;
             this.Pts.P = P;
        end      
    end
        
    
    methods (Access = public)
        
        function Int =  ComputeIntegrationVector(this)
            [~,wInt] = LGLPointsWeights(this.N-1);
            [~,J] = PhysSpace(this,this.Pts.x);
            
            J(J==inf)  = 0;  % assume integration well-defined
            J(J==-inf)  = 0;  % assume integration well-defined
            J(isnan(J))= 0;  % at infinity
            
            Int = wInt.*(J.');
            this.Int = Int;
        end    
        function Diff = ComputeDifferentiationMatrix(this)
        
            %CompDiff  = barychebdiff(this.Pts.x,4);   
            CompDiff  = LGLDiff(this.Pts.x,this.Pts.P,2);   
            Diff = PhysicalDerivatives_1D(@this.PhysSpace,this.Pts.x,CompDiff);            
          %  [~,~,dx,ddx] = PhysSpace(this,this.Pts.x); 
%            Diff.Dy   = diag(dx)*CompDiff.Dx;  
%            Diff.DDy  = diag(ddx)*CompDiff.Dx + (diag(dx.^2))*CompDiff.DDx;    
            this.Diff = Diff;
        
        end        
        function Interp = ComputeInterpolationMatrix(this,interp,saveBool)
                       
            InterPol = barychebevalMatrix(this.Pts.x,interp);
            Interp = struct('InterPol',InterPol,...
                            'Nplot',length(interp),...
                            'N',this.N,...
                            'pts',PhysSpace(this,interp),...
                            'xPlot',interp,...
                            'yPlot',this.PhysSpace(interp));

            if((nargin >= 3) && saveBool)
                this.Interp = Interp;
            end  
        end      
        function Ind = ComputeIndices(this)
            Ind = GetIndicesInterval(this);
            this.Ind = Ind;
        end       
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            N  = this.N; 
            Pts = this.Pts;
            Int = this.Int;  % 1 x N

            % find size of function matrix by making a dummy version
            fPTemp = f(Pts.y-Pts.y);
            fDim   = size(fPTemp);
            nElts  = prod(fDim(2:end)); % first dimension stores values
            
            IntT = Int.';  % N x 1
            
            IntT = IntT(:,ones(1,nElts)); % N x nElts
            IntT = reshape(IntT,fDim);    % size(f)
            
            M_conv = zeros([N,N,fDim(2:end)]);
            
            Mmask = repmat({':'},[1,fDim]);
            
            for i=1:N 
                fP          = f(Pts.y(i)-Pts.y);  
                Mmask{1} = i;
                M_conv(Mmask{:}) = IntT.*fP;
            end
            M_conv(isnan(M_conv)) = 0;
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
            
        end                
    end
    
end