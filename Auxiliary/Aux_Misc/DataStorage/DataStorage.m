function [Data,recompute,Parameters] = DataStorage(nameDir,func,Parameters,OtherInputs,recompute,ignoreList,addLoadTime,skipDataStorage)
%Input:
%Name of Folder,List of Parameters
%if file is to be written: Data
%
%!!! Attention: As Parameters, numbers will be trimmed to 12 digits !!!

if(nargin == 8 && skipDataStorage)
    if (nargin(func) >= 2)
        Data = func(Parameters, OtherInputs);
    else
        % func may be using name-value arguments rather than structs
        Data = CallNameValueFunc(func, Parameters, OtherInputs);
    end
    recompute = true;
    return
end

%%*****************************
%Check if there is a personal user!
global recomputeAll loadAll
if(isfield(OtherInputs,'recomputeAll'))
    recomputeAll =  OtherInputs.recomputeAll;
end
if(isfield(OtherInputs,'loadAll'))
    loadAll =  OtherInputs.loadAll;
end

%*****************************
%Trim numbers of parameters
OriginalParameters = Parameters;

%Parameters         = FlattenStructure(Parameters,10,'AllStr');
Parameters         = FlattenStructure(FlattenStructure(Parameters),10,'AllStr');

[h1s,Parameters.Function] = fileparts(func2str(func));
if(isfield(OriginalParameters,'Identifier'))
    identifier = OriginalParameters.Identifier;
else
    identifier = [];
end
if(isfield(Parameters,'Comments'))
    comments            = Parameters.Comments;
else
    comments = '';
end
hname               = getFilename(Parameters);
filename            = [hname,'.mat'];
fileParamTxtname    = [hname,'.txt'];
%*****************************

global dirData
if(isfield(OtherInputs,'dirData'))
    dirDataTemp = OtherInputs.dirData;
else
    dirDataTemp = dirData;
end
old_dirData = dirDataTemp;

if(~isempty(nameDir))
    DataFolder = [dirDataTemp filesep nameDir];
else
    DataFolder = dirDataTemp;
end

index      = LoadIndexFiles(DataFolder);

if(nargin <= 4)
    recompute = [];
end

if(nargin < 6)
    ignoreList = {};
end

if(nargin < 7)
    addLoadTime = true;
end

%1st Step: Search for File and load if found
if((nargin >= 5) && ~(isempty(recompute)) &&  ischar(recompute))
    load(recompute);
    disp(['Stored data from ',recompute,' will be used..']);
    recompute = false;
elseif((nargin >= 5) && ~(isempty(recompute)) && (recompute == 2))
    %(1)First Option: Select Input File
    [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);
    load([DataFolder,FileName]);
    disp(['Stored data from ',[DataFolder,FileName],' will be used..']);
    recompute = false;
else
    %(2) Second Option: Search file
    
    i          = DataStorage_FindFile(index,Parameters,ignoreList);
    
    if(~isempty(i))
        Parameters = index{i};
        
        [~,filename,ext]    = fileparts(index{i}.Filename);
        filename            = [filename,ext];
        % filename            = index{i}.Filename;
        Parameters.Filename = [DataFolder filesep filename(1:end-4)];
        
        fileParamTxtname    = [filename(1:end-4),'.txt'];
        
        try
            load([DataFolder filesep filename]);
            disp(['Data file found in ',DataFolder filesep filename,' ...']);
            
            if(isfield(OtherInputs,'recomputeAll'))  % mostly for parfor
                recompute = OtherInputs.recomputeAll;
            elseif(recomputeAll)
                recompute = true;
            elseif(loadAll)
                recompute = false;
            elseif((isempty(recompute))) %ask if value of recompute is not given as an input
                no = fprintf('Do you want to load data (press "l"), recompute Data (press any other key), or wait for 2 seconds.');
                if(getkeywait(2) == 108)
                    [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);
                    load([DataFolder,FileName]);
                    disp(['Stored data from ',[DataFolder,FileName],' will be used..']);
                    
                    recompute = false;
                    
                elseif(getkeywait(2) == -1)
                    for ih = 1:no
                        fprintf('\b');
                    end
                    recompute = false;
                else
                    for ih = 1:no
                        fprintf('\b');
                    end
                    no = fprintf('Thanks. Data will be recomputed.\n');
                    recompute = true;
                end
            end
            if(~recompute)
                if(addLoadTime)
                    %append configuration to index file
                    fileID = fopen([DataFolder filesep fileParamTxtname],'a');
                    fprintf(fileID,'# %s \n',['loaded on ' datestr(clock) ' : ' comments]);
                    fclose(fileID);
                end
            end
        catch mexception
            if mexception.identifier == "MATLAB:load:couldNotReadFile"
                disp('No mat file corresponding to found text file. Recomputing data...');
                recompute = true;
            else
                rethrow(mexception);
            end
        end
    end
    
    if(~exist('Data','var'))
        if(isfield(OtherInputs,'SkipLoadQuestion') && ~OtherInputs.SkipLoadQuestion)
            no = fprintf('Do you want to load data (press "l"), or wait for 2 seconds.');
            if(getkeywait(2) == 108)
                [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);
                load([DataFolder,FileName]);
                disp(['Stored data from ',[DataFolder,FileName],' will be used..']);
                
                Parameters.Filename =  [DataFolder,FileName];
                
                recompute = false;
            else
                for ih = 1:no
                    fprintf('\b');
                end
                recompute = true;
            end
        else
            recompute = true;
        end
    end
    
    %2b) If File not found: Compute Data and Save Data
    %if(~exist('Data','var') || ((nargin == 5) && recompute == true))
    if(recompute)
        startComp = tic;
        fprintf(['Computing data, starting at ', datestr(now), '...\n']);
        if (nargin(func) >= 2)
            if (nargout(func) >= 2)
                [Data,Parameters.Results] = ...
                    func(OriginalParameters,OtherInputs);
            else
                Data = func(OriginalParameters, OtherInputs);
            end
        else
            % func may be using name-value arguments rather than structs
            if (nargout(func) >= 2)
                error('Not implemented yet');
            else
                Data = CallNameValueFunc(func, OriginalParameters, OtherInputs);
            end
        end
        t                               = toc(startComp);
        
        if(~exist(DataFolder,'dir'))
            disp('Folder not found. Creating new path..');
            mkdir(DataFolder);
        end
        
        
        save([DataFolder filesep filename],'Data','-v7.3');
        
        %            lastwarn('');
        %
        %            save([DataFolder filesep filename],'Data');
        %
        %            % deal with case when file is >2GB
        %            [~,wid]=lastwarn;
        %            if(strcmp(wid,'MATLAB:save:sizeTooBigForMATFile'))
        %               save([DataFolder filesep filename],'Data','-v7.3');
        %            end
        
        disp(['Data saved in ',[DataFolder filesep filename],'.']);
        
        Parameters.Results.comptime_sec = t;
        Parameters.Results.comp_at      = datestr(now);
        if(isnumeric(Data) && isscalar(Data) && (Data == 0))
            disp(['Error after ',sec2hms(t), ' (hrs:min:sec) ']);
            comments = ['Computation failed! ',comments];
        else
            disp(['Data recomputed: ',sec2hms(t),' (hrs:min:sec)']);
        end
        Parameters.Results.comments  = comments;
        Parameters.Filename          = filename;
        
        Struct2File([DataFolder filesep fileParamTxtname],Parameters,...
            ['Computed at: ',datestr(now),...
            ' Computation time: ',sec2hms(t), ' (hrs:min:sec) ',comments]);
        
        %disp(['Index file saved in ',[DataFolder filesep fileParamTxtname],'..']);
        Parameters.Filename = [DataFolder filesep Parameters.Filename(1:end-4)];
    end
end
ChangeDirData(old_dirData);


    function name = getFilename(params)
        if(isempty(identifier))
            name   = '';
        elseif(ischar(identifier))
            identifier = {identifier};
        end
        if(iscell(identifier))
            name = '';
            nIdent = length(identifier);
            for iIdent = 1:nIdent
                if(isfloat(params.(identifier{iIdent})))
                    name = [name identifier{iIdent},'_',num2str(params.(identifier{iIdent})),'_'];
                else
                    name = [name identifier{iIdent},'_',params.(identifier{iIdent}),'_'];
                end
            end
        end
        name   = ([name,getTimeStr()]);
    end

%     function name = getFilename(params)
%         if(isempty(identifier))
%             name   = '';
%         elseif(isfloat(params.(identifier)))
%             name = [identifier,'_',num2str(params.(identifier)),'_'];
%         else
%             name = [identifier,'_',params.(identifier),'_'];
%         end
%         name   = ([name,getTimeStr()]);
%     end

end

% TODO: Make this work for more than one output argument
function Data = CallNameValueFunc(func, Parameters, OtherInputs)
% Call a function with name-value pair arguments
% Convert struct to cell array
P = Parameters;
OI = OtherInputs;
% Strip off DataStorage specific fields
for field = ["Function", "Filename", "Results", "Comments", "Identifier"]
    if isfield(P, field)
        P = rmfield(P, field);
    end
end
for field = ["recomputeAll", "loadAll", "dirData", "SkipLoadQuestion"]
    if isfield(OI, field)
        OI = rmfield(OI, field);
    end
end

Pcell = namedargs2cell(P);
OIcell = namedargs2cell(OI);
% Call the function
Data = func(Pcell{:}, OIcell{:});
end
