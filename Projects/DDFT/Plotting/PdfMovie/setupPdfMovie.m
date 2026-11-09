function [nd,pdfFileNames] = setupPdfMovie(outTimes,pdfDir)

    nDigits=ceil(log10(length(outTimes)));
    nd=['%0' num2str(nDigits) 'd'];
    
    pdfFileNames = [];
    
    if ~exist(pdfDir,'dir')
        mkdir(pdfDir);
    end
end
