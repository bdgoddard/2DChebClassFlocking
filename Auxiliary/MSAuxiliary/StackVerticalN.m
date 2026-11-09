function s = StackVerticalN(varargin)
    s = [];
    for i = 1:length(varargin)
        si = varargin{i};
        s = [s;si];
    end

end