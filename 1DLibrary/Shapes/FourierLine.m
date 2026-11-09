classdef FourierLine < Fourier
    
    properties
        yMin
        yMax
        L
    end
    
    methods
        function this = FourierLine(Geometry)
            this@Fourier(Geometry.N);
            if(isfield(Geometry,'yMin'))
                this.yMin = Geometry.yMin;
            else
                this.yMin =0;
            end
            if(isfield(Geometry,'yMax'))
                this.yMax = Geometry.yMax;
            else
                this.yMax =2*pi;
            end

            this.L = this.yMax - this.yMin;
            this.polar = 'cart';
                                    
            InitializationPts(this);  
        end        
    end
    
   methods (Access = public)      
        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace(this,x)
            [y,dy,dx,ddx,dddx,ddddx] = LinearMap01(x,this.yMin,this.yMax);
        end             
        function xf = CompSpace(this,y)
            xf = InvLinearMap01(y,this.yMin,this.yMax);
        end        
                
        function M_conv = ComputeConvolutionMatrix(this,f,shapeParams,parent)
            %{
            Strategy:
            We wish to compute \int_D(x) f(x-y) g(y) dy, for some domain D
            (e.g. disc, annulus) centred at x.
            Let z = y-x then dy = dz, y = x + z and 
            y \in D(x) => z + x \in D(x) or z \in D(0) and we get
            \int_D(0) f(-z) g(x+z) dz
            
            In matlab this is
            M_conv = Int_D(0) .* f(-z) .* IP_D(x)
            %}            
            
            if((nargin == 4) && parent)
                M_conv = ComputeConvolutionMatrix@Fourier(this,f);
                return;
            end

            if(nargin<3)
                shapeParams = struct;
            end
                        
            % disp('Computing Convolution matrices...'); 
            
            if(isfield(shapeParams,'yMin') && isfield(shapeParams,'yMax'))
                % subdomain
                nSubShapes = 1;
                
                geom.N    = shapeParams.N;
                geom.yMin = shapeParams.yMin;
                geom.yMax = shapeParams.yMax;
                subShape  = SpectralLine(geom);
              
            elseif(isfield(shapeParams,'yMax'))
                % subdomain
                nSubShapes = 1;
                
                geom.N    = shapeParams.N;
                geom.yMin = -shapeParams.yMax;
                geom.yMax = shapeParams.yMax;
                subShape  = SpectralLine(geom);
                
            else
                                
                %whole interval
                M_conv = ComputeConvolutionMatrix@Fourier(this,f);
                return;
                
            end
            
            % initialize matrix
            fP = f(-subShape(1).Pts.y);

            nDim = length(size(fP));
            if(nDim>2)
                fP_old = fP;
                fP = reshape(fP,size(fP,1),[]);
            end

            nf = size(fP,2);
            M_conv = zeros(this.N,this.N,nf);  
           
            for iSubShape = 1:nSubShapes
             
                ySub = subShape(iSubShape).Pts.y;
             
                % Integration and interpolation on D(0)
                Int_i    = subShape(iSubShape).ComputeIntegrationVector();
                
                % f(-z)
                fP = f(-ySub);
                
                if(nDim>2)
                    fP_old = fP;
                    fP = reshape(fP,size(fP,1),[]);
                end
                IntRep = Int_i(ones(1,nf),:);

                for i=1:this.N
                
                    Y0 = this.Pts.y(i); 
                    % Compute points in D(x) = D(y0))
                    subConvShapePts.y = subShape(iSubShape).Pts.y+Y0;
                    
                    % and the interpolation onto D(x) from the full domain
                    IP          = SubShapePts(this,subConvShapePts);
                                   
                    IntFTemp = (IntRep.*fP.')*IP;

                    if(size(IntFTemp,1)>1)
                        IntFTemp = IntFTemp.';
                    end
                    
                    M_conv(i,:,:) = squeeze(M_conv(i,:,:)) + IntFTemp;

                end

                %close(hw);
                
                M_conv(isnan(M_conv)) = 0;
                            
                if(nDim>2)
                    nfTemp = size(fP_old);
                    nfTemp = nfTemp(2:end);
                    M_conv = reshape(M_conv,[ this.N, this.N, nfTemp ]);
                end

            end % subshape loop
                
            % disp('... done.'); 
            
        end % convolution               
        
        function yInt = MapToInterval(this,y)
            yInt = mod(y,this.L)+this.yMin;
        end
        
   end
    
end