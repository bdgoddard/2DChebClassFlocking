classdef (Abstract) FiniteDifferenceFiniteDifference < Shape
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = FiniteDifferenceFiniteDifference(N1,N2)
             this@Shape(N1,N2);                                    
             
             this.Pts.x1 = FDPointsWeights(N1);
             this.Pts.x2 = FDPointsWeights(N2);
             
             this.Pts.N1 = N1;
             this.Pts.N2 = N2;
        end      
    end
    
    %**********************************************
    %************   Initializations    *************
    %**********************************************
    methods (Access = public) 
        
        function int = ComputeIntegrationVector(this)
            [~,dy1]        = this.PhysSpace1(this.Pts.x1);
            [~,dy2]        = this.PhysSpace2(this.Pts.x2);             
            
            [~,wInt1] = FDPointsWeights(this.N1);
            [~,wInt2] = FDPointsWeights(this.N2);
            
            this.Int = kron( dy1'.*wInt1, dy2'.*wInt2);            
            int = this.Int;
        end   
        
        function Diff = ComputeDifferentiationMatrix(this,pts)
            
            if(nargin == 1)
                pts = this.Pts;                
            end
                        
            Diff1 = FDDiff(this.N1);    
            Diff2 = FDDiff(this.N2);
            
            Sel = {'Dy1' ;'DDy1' ; 'Dy2'; 'DDy2' ;...
                  'Dy1Dy2'; 'Lap' ;'grad' ;'div';...
                  'gradDiv'; 'LapVec'}; 

            Diff = PhysicalDerivatives(this,pts,Sel,Diff1,Diff2,2);        

        end              
        function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)
            Interp1     = FDInterpolation(this.Pts.x1,interp1); 
            Interp2     = FDInterpolation(this.Pts.x2,interp2); 


            InterPol  = kron(Interp1,Interp2);
            
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
            Ind      = GetIndicesBox(this);
            this.Ind = Ind;
        end    
        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)

            N1 = this.N1;   N2 = this.N2;

            xs  = this.Pts.y1_kv; 
            ys  = this.Pts.y2_kv;            
            M_conv = zeros(N1*N2,N1*N2);

            for i=1:(N1*N2)
                M_conv(i,:)  = this.Int.*f(xs(i)-xs,ys(i)-ys)';
            end   
            if((nargin== 3) && saveBool)
                this.Conv = M_conv;
            end
            
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


