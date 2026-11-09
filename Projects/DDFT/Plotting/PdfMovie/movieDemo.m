function movieDemo()

%--------------------------------------------------------------------------
% Choose what type of movie to make

doPdfs = false; % pdf movie for beamer
doMp4 = true;   % mp4 movie for beamer
%--------------------------------------------------------------------------

% set up a spectral line for the plotting
geom.N = 20;
geom.yMin = 0;
geom.yMax = 1;

aLine = SpectralLine(geom);
aLine.ComputeAll;
aLine.ComputeInterpolationMatrix([-1:0.01:1]',true);

y = aLine.Pts.y;

% plot time
tMax = 1;
plotTimes = 0:tMax/100:tMax;
nPlots = length(plotTimes);

% options for the two plots
optsPlot1.plain = true; % removes green circles at collocation points
optsPlot1.linecolor = 'r';
optsPlot1.linestyle = '-';
optsPlot1.linewidth = 2;

optsPlot2.plain = true;
optsPlot2.linecolor = 'b';
optsPlot2.linestyle = '--';
optsPlot2.linewidth = 1;

% set up figure and axes
hf = figure('Position',[100,100,1200,800]);
ha = axes;

% create a directory to hold the pdf files
if(doPdfs)
    pdfDir = [pwd filesep 'TempMovies' filesep];
    [nd,pdfFileNames] = setupPdfMovie(plotTimes,pdfDir);
end

% open an mp4 file to write to
if(doMp4)
    v = VideoWriter('movieDemo.mp4','MPEG-4');
    v.FrameRate=10;
    v.Quality = 100;
    open(v)
end

% do the plotting
for iPlot = 1:nPlots
    tPlot = plotTimes(iPlot);
    f1 = toPlot1(y,tPlot);
    f2 = toPlot2(y,tPlot);
    
    %plot with corresponding options
    aLine.plot(f1,optsPlot1);
    hold on
    aLine.plot(f2,optsPlot2);
    hold off
    
    % plot setup
    title(['Time = ' num2str(tPlot)]);
    set(gca,'ylim',[-1,1]);
    
    % save the pdf in a file with a name corresponding to the frame number
    if(doPdfs)
        pdfFileNames = savePdfFrame(pdfDir,pdfFileNames,iPlot,nd,hf);
    end
    
    % save the mp4 frame
    if(doMp4)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    pause(0.1)
    
end

% combine pdf files into a single one
if(doPdfs)
    makePdfMovie(pdfDir,pdfFileNames);
end

% close & write the mp4 file
if(doMp4)
    close(v);
end

% simple functions for plotting

    function f = toPlot1(y,t)
        f = exp(-(y-t).^2/0.1);
    end

    function f = toPlot2(y,t)
        f = t*sin(2*pi*y);
    end


end