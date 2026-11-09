function Conv = ComputeConvolutionWithParameters(domain,opts)

    convFn = @convolutionWrapper;
    recompute = false;
    optsOther.SkipLoadQuestion = true;

    convDir = opts.convDir;
    convParams.geom = opts.geom;
    convParams.optsKernel = opts.optsKernel;

    if(isfield(opts,'convFn'))
        convParams.convFn = opts.convFn;
        if(isfield(opts,'optsConvFn'))
            convParams.optsConvFn = opts.optsConvFn;
        end
    end
    
    dataConv = DataStorage(convDir,convFn,convParams,optsOther,[],[],recompute);
    
    Conv = dataConv.Conv;

    
    function data = convolutionWrapper(params,otherInputs)
        
        if(isfield(params,'convFn'))
            
            if(isfield(params,'optsConvFn'))
                data.Conv = domain.(params.convFn)(@kernelWrapper,params.optsConvFn);
            else
                data.Conv = domain.(params.convFn)(@kernelWrapper);
            end
        else
            data.Conv = domain.ComputeConvolutionMatrix(@kernelWrapper);
        end

        function K = kernelWrapper(y)
            kernelFn = str2func(params.optsKernel.kernelFn);
            if(nargin(kernelFn) == 1)
                K = kernelFn(y);
            else
                K = kernelFn(y,params.optsKernel);
            end
        end

    end

end