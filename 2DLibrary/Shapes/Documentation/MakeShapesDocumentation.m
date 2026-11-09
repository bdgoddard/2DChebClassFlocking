function MakeShapesDocumentation(opts)
    % Make documentation for Shapes in 2DChebClass.  Produces an output
    % with an example geometry plot and the help text for each class that
    % supports this.
    %
    % Inputs:
    %  opts - a structure containing (all optional)
    %    - outputDir - the directory to save the file in
    %    - outputFormat - one of 'html', 'doc', 'latex', 'ppt', 'xml',
    %                     'pdf' (default)
    %    - alphabetical - whether to sort by file name

    if(nargin<1)
        opts = struct;
    end
   
    % set output dir
    if(isfield(opts,'outputDir'))
        outputDir = opts.outputDir;
    else
        outputDir = ['.' filesep 'Library' filesep '2DLibrary' filesep 'Shapes' filesep 'Documentation' filesep 'Output'];
    end
    
    % set format
    if(isfield(opts,'format'))
        outputFormat = opts.outputFormat;
    else
        outputFormat = 'pdf';
    end

    if(isfield(opts,'alphabetical'))
        alphabetical = opts.alphabetical;
    else
        alphabetical = false;
    end

    codeToEvaluate = ['ShapesDocumentation(' num2str(alphabetical) ');'];
    publish('ShapesDocumentation','codeToEvaluate',codeToEvaluate, ...
        'showCode',false,'format',outputFormat,'outputDir',outputDir);

    
end