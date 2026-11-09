function makePdfMovie(pdfDir,pdfFileNames)

    fprintf(1,'Combining pdf ... ');
    fullPdfFile=[pdfDir  'movie.pdf'];

    gs = getGS();

    gsCmd= [gs ' -dNOPAUSE -sDEVICE=pdfwrite ' ...
              '-sOUTPUTFILE=' fullPdfFile ' -dBATCH -dQUIET ' pdfFileNames];

    system(gsCmd);
    fprintf(1,'Finished\n');
    
end