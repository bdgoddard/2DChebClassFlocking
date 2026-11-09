function [bincounts,edges] = PlotKuramotoSDE(data,optsPlot)

if(nargin==1)
    optsPlot = struct;
end

t = data.tStoc;
x = data.x;
opts = data.opts;

nParticles = opts.nParticles;

if(length(size(x))==3) % there are multiple runs
    nRuns = size(x,3);
else
    nRuns = size(x,1)/nParticles; % if sampling
end

if(isfield(optsPlot,'doPlots'))
    doPlots = optsPlot.doPlots;
else
    doPlots = true;
end

if(doPlots)
    if(isfield(optsPlot,'ha'))
        ha = optsPlot.ha;
    else
        figure('Position',[100,100,800,800])
        ha = axes;
    end
end


if(isfield(optsPlot,'nBins'))
    nBins = optsPlot.nBins;
else
    nBins = round(nParticles/10);
end

if(isfield(optsPlot,'xLim'))
    xLim = optsPlot.xLim;
else
    xLim = [0,2*pi];
end

dx = (xLim(2)-xLim(1))/nBins;
edges = xLim(1):dx:xLim(2);

if(isfield(optsPlot,'tPlot'))
    [~,tPos] = min(abs(t-optsPlot.tPlot));
else
    % default to final time
    tPos = length(t);
end

xHist = x(:,tPos,:);
xHist = xHist(:);

[bincounts,edges] = histcounts(xHist,edges);

bincounts = bincounts/nParticles/nRuns/dx;

if(isfield(optsPlot,'linecolor'))
    linecolor = optsPlot.linecolor;
else
    linecolor = 'r';
end

if(doPlots)
    hs = scatter(ha,edges(1:end-1)+dx/2,bincounts);
    set(hs,'MarkerEdgeColor',linecolor);
    set(hs,'LineWidth',1.5);
    if(isfield(optsPlot,'handlevisibility'))
        set(hs,'HandleVisibility',optsPlot.handlevisibility);
    end
    
    xlim(ha,xLim);

    if(isfield(optsPlot,'yLim'))
        ylim(ha,optsPlot.yLim);
    end

    if(isfield(optsPlot,'xLabel'))
        xlabel(ha,optsPlot.xLabel,'interpreter','latex');
    else
        xlabel(ha,'$\theta$','interpreter','latex');
    end
    if(isfield(optsPlot,'yLabel'))
        ylabel(ha,optsPlot.yLabel,'interpreter','latex');
    else
        ylabel(ha,'Density','interpreter','latex');
    end
    
end

end