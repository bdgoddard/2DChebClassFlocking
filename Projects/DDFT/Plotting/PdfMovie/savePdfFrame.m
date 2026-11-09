function pdfFileNames = savePdfFrame(pdfDir,pdfFileNames,iTimes,nd,hf)

    outputFile=[pdfDir num2str(iTimes,nd) '.pdf'];
    % add to list of pdf files used to make movie
    pdfFileNames = cat(2, pdfFileNames, [' ' outputFile]);
    save2pdf(outputFile,hf,[]);
end
