function [hf,ha,plotFile] = DocumentationPlotDefault(className,opts)
    % set a sensible default geometry and make the documentation
    % plot
    
    % set default line width and text sizes
    set(groot,'defaultlinelinewidth', 2);
    set(groot,'defaulttextfontsize', 20);
    set(groot,'defaultaxesfontsize', 20);

    % deal with empty options
    if(nargin < 2)
        opts = struct;
    end
    
    % get the default geometry by running the static method
    try
        geom = eval([className '.DocumentationPlotDefaultGeometry']);
        %classFnTest = str2func([className '.DocumentationPlotDefaultGeometry']);
    catch
        hf = [];
        ha = [];
        plotFile = '';
        return
    end

    % make an instance of the class using the geometry
    classFn = str2func(className);
    if(nargin(classFn) == 0) % specified by argument list
        aClassObject = classFn(geom{:});
    else % specified by structure
        aClassObject = classFn(geom);
    end

    % make the plot
    [hf,ha,plotFile] = aClassObject.DocumentationPlot(opts);

end
