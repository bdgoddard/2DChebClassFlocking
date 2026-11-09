pdfFileNames=[];
pdfDir = ['TempMovie' filesep];

if(~exist(pdfDir,'dir'))
    mkdir(pdfDir);
end
    
nPlots = 10;

nDigits=ceil(log10(nPlots));
nd=['%0' num2str(nDigits) 'd'];

x = 0:0.01:1;

figure;

for iPlot = 1:nPlots

    plot(x,x.^iPlot);

    outputFile=[pdfDir num2str(iPlot,nd) '.pdf'];
    % add to list of pdf files used to make movie
    pdfFileNames = cat(2, pdfFileNames, [' ' outputFile]);

    % save the pdf file for this time step
    save2pdf(outputFile,gcf);
end

fprintf(1,'Combining pdf ... ');
fullPdfFile=[pdfDir  'movie.pdf'];

switch computer
		case {'MAC','MACI','MACI64'}			
            gs= '/usr/local/bin/gs';
		case {'PCWIN','PCWIN64'}
            gs= 'gswin32c.exe';
        otherwise
            gs= 'gs';
end

gsCmd= [gs ' -dNOPAUSE -sDEVICE=pdfwrite ' ...
          '-sOUTPUTFILE=' fullPdfFile ' -dBATCH -dQUIET ' pdfFileNames];

system(gsCmd);
fprintf(1,'Finished\n');