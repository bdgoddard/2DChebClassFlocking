function [dataPDE,dataSDE] = KuramotoComparison(opts,optsPlot)

% Deal with any missing inputs
if(nargin<1)
    opts = struct;
end
if(nargin<2)
    optsPlot = struct;
    optsPlot.doPlots = false;
end

% Determine which computations we're doing
if(isfield(opts,'doPDE'))
    doPDE = opts.doPDE;
else
    doPDE = true;
end

if(isfield(opts,'doSDE'))
    doSDE = opts.doSDE;
else
    doSDE = true;
end

% Check if we are under or overdamped
if(isfield(opts,'phasespace'))
    phasespace = opts.phasespace;
else
    % default is overdamped
    phasespace = false;
    opts.phasespace = phasespace;
end

% Do the relevant PDE dynamics
if(doPDE)
    if(phasespace)
        dataPDE = KuramotoPhaseSpace(opts,optsPlot);
    else
        dataPDE = KuramotoOverdamped(opts,optsPlot);
    end
else
    dataPDE = [];
end

% Do ther relevant SDE dynamics
if(doSDE)
    % automatically switches between under- and over-damped
    dataSDE = KuramotoSDEDynamics(opts,optsPlot);
else
    dataSDE = [];
end

% Make a movie
if(optsPlot.doMp4)
    KuramotoMovie(dataPDE,dataSDE,optsPlot);
end

% Do snapshot plot
if(optsPlot.doSnapshots)
    KuramotoSnapshots(dataPDE,dataSDE,optsPlot);
end


end