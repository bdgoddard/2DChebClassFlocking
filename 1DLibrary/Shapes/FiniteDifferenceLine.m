classdef FiniteDifferenceLine < FiniteDifference
    
    properties
        yMin
        yMax
    end
    
    methods
        function this = FiniteDifferenceLine(Geometry)
            this@FiniteDifference(Geometry.N);
            this.yMin = Geometry.yMin;
            this.yMax = Geometry.yMax;
            
            this.polar = 'cart';
            
            InitializationPts(this);  
        end
        
    end
    
    
   methods (Access = public)
       
        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace(this,x)
            [y,dy,dx,ddx,dddx,ddddx] = LinearMap(x,this.yMin,this.yMax);
        end
        
        function xf = CompSpace(this,y)
            xf  = InvLinearMap(y,this.yMin,this.yMax);
        end                        
        
        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams,parent)
            error('Finite difference convolution not implemented yet');
        end
        
    end
    
end