function Struct2File(filename,st,comment_str)
%Saves elements of Structure in File
%Numbers will be saved up to 14 digits

if((nargin < 3) && isfield(st,'SecondLine'))
    comment_str = st.SecondLine;
    st          = rmfield(st,'SecondLine');
end

st = FlattenStructure(st);


fileID = fopen(filename,'w');

%Print comment:
fprintf(fileID,'# ************************************************************ \n');
fprintf(fileID,'# %s \n',comment_str);
fprintf(fileID,'# ************************************************************ \n');


names = fieldnames(st);

for i = 1:length(names)
    str = st.(names{i});
    if isnumeric(str)
        str = num2str(str,12);
    elseif iscell(str)
        str = cell2str(str);
    elseif isa(str, 'function_handle')
        str = char(str);
    elseif isstring(str)
        str = sprintf('%s ', str);
    end
    fprintf(fileID,'%s : %s\n',names{i},str);
end
fclose(fileID);

end