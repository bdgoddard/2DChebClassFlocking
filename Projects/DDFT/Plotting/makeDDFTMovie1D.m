function makeDDFTMovie1D(DDFT)

    shape = DDFT.shape;
    rho = DDFT.rho_t;
    flux = DDFT.flux_t;
    outTimes = DDFT.outTimes;
    
    nPlots = length(outTimes);
    
    movieDir = ['Results' filesep 'Movies' filesep];
    
    if(~exist(movieDir,'dir'))
        mkdir(movieDir);
    end
    
    pdfFileNames = [];
    nDigits=ceil(log10(length(outTimes)));
    nd=['%0' num2str(nDigits) 'd'];    
    
    for iPlot = 1:nPlots
        rhot = rho(:,:,iPlot);
        subplot(2,1,1)
        hold off
        shape.plot(rhot)
        ylabel('$\rho$', 'Interpreter','LaTeX');
        xlabel('$x$', 'Interpreter','LaTeX');
        title(['t = ' num2str(outTimes(iPlot))]);
        subplot(2,1,2)
        hold off
        fluxt = flux(:,:,iPlot);
        shape.plot(fluxt)
        ylabel('Flux', 'Interpreter','LaTeX');
        xlabel('$x$', 'Interpreter','LaTeX');
        outputFile=[movieDir num2str(iPlot,nd) '.pdf'];
        % add to list of pdf files used to make movie
        pdfFileNames = cat(2, pdfFileNames, [' ' outputFile]);
        
        save2pdf(outputFile,gcf);
        
    end

    
    switch computer
		case {'MAC','MACI','MACI64'}			
            gs= '/usr/local/bin/gs';
		case {'PCWIN','PCWIN64'}
            gs= 'gswin32c.exe';
        otherwise
            gs= 'gs';
    end

    fprintf(1,'Combining pdf ... ');
    fullPdfFile=[movieDir  'movie.pdf'];

    gsCmd= [gs ' -dNOPAUSE -sDEVICE=pdfwrite ' ...
              '-sOUTPUTFILE=' fullPdfFile ' -dBATCH -dQUIET ' pdfFileNames];

    system(gsCmd);
    fprintf(1,'Finished\n');
    
end