function Conv = ConvMatrices(V2params,MS)
   
    V2 = str2func(V2params.V2Str);
    % Computing Convolution Matrix    
    Conv = MS.ComputeConvolutionMatrix(@KernelUniversal,true);  

%--------------------------------------------------------------------------
% The below determines whether the V2 function takes parameters or not.
% This needs to be done to be able to datastore the convolution
% computations...

	function z = KernelUniversal(r)                 
        if(nargin(V2)==1)
            z = V2(r);   
        else
            z = V2(r,V2params);   
        end
    end

end

% This function allows to compute convolutions and pass parameters to
% the function V2
% 'based on' FexMatrices_Meanfield.m