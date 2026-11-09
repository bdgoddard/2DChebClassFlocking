function gs = getGS()

switch computer
		case {'MAC','MACI','MACI64','MACA64'}			
            gs= '/usr/local/bin/gs';
		case {'PCWIN','PCWIN64'}
            gs= 'gswin32c.exe';
        otherwise
            gs= 'gs';
end

end
