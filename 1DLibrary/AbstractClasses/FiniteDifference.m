classdef (Abstract) FiniteDifference < Interval
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = FiniteDifference(N)
             this@Interval(N);             
             this.Pts.x = FDPointsWeights(N); 
        end      
    end
        
    
    methods (Access = public)
        
        function Int =  ComputeIntegrationVector(this)
            [~,wInt] = FDPointsWeights(this.N);
            [~,J] = PhysSpace(this,this.Pts.x);
                        
            Int = wInt.*(J.');
            this.Int = Int;
        end 
        
        function Diff = ComputeDifferentiationMatrix(this)
        
            CompDiff  = FDDiff(this.N);   
            Diff = PhysicalDerivatives_1D(@this.PhysSpace,this.Pts.x,CompDiff);            
            this.Diff = Diff;
        
        end
        
        function Interp = ComputeInterpolationMatrix(this,interp,saveBool)
            
            InterPol = FDInterpolation(this.Pts.x,interp);
            Interp = struct('InterPol',InterPol,...
                            'Nplot',length(interp),...
                            'N',this.N,...
                            'pts',PhysSpace(this,interp));

            if((nargin >= 3) && saveBool)
                this.Interp = Interp;
            end  
        end      
        
        function Ind = ComputeIndices(this)
            Ind = GetIndicesInterval(this);
            this.Ind = Ind;
        end       
        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)

            error('Finite difference convolution not implemented yet');

        end
    end
    
end