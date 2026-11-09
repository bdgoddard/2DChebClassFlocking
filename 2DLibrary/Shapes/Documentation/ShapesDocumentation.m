function ShapesDocumentation(alphabetical)
    % Make documentation for Shapes in 2DChebClass and displays on screen.
    % Can be used with MakeShapesDocumentation to produce file-based
    % output.
    %
    % Inputs:
    %   - alphabetical - whether to sort by file name

    if(nargin<1)
        alphabetical = false;
    end
    
    fileNames = {};  % list of .m files
    classNames = {}; % list of clas names
    md = [];  % metadata to determine which .m files are classes
    
    % directory with all of the Shapes
    ShapesDir = ['.' filesep 'Library' filesep '2DLibrary' filesep 'Shapes' filesep];
    
    % get all m-files in any subdirectories
    listingM = dir(fullfile(ShapesDir, ['**' filesep '*.m'])); 
    nFiles = length(listingM);
    
    % get a list of all file names without path or extension
    for iFile = 1:nFiles
        fileName = listingM(iFile).name(1:end-2);
        fileNames = [fileNames, fileName];  %#ok
    end
    
    % establish which are valid classes
    for iFile = 1:nFiles
        fileName = fileNames{iFile};
        metaDataCmd = ['md  = ?' fileName ';'];
        eval(metaDataCmd); % is empty if .m file isn't a class
    
        if(~isempty(md))
            classNames = [classNames, fileName];  %#ok
        end
    end
    
    % put into alphabetical order
    if(alphabetical)
        [~,sortOrder] = sort(lower(classNames));
        classNames = classNames(sortOrder);
    end
    nClasses = length(classNames); % number of classes, <= number of files
    
    % produce figure and help text
    for iClass = 1:nClasses
    
        % do the default documentation plot
        className = classNames{iClass};
        hf = DocumentationPlotDefault(className);
    
        % only display anything if there's a default plot
        if(~isempty(hf))
            snapnow;  % take snapshot if default plot
            close(hf); 
    
            % get help text and strip off standard ending
            helpText = help(className);
            helpText = eraseBetween(helpText,'Documentation for',length(helpText),'Boundaries','inclusive');
            
            % print the class name and help text
            disp(className);
            disp(helpText);
        end
    end

end