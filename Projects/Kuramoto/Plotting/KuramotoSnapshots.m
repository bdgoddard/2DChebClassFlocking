function KuramotoSnapshots(dataPDE,dataSDE,optsPlot)

    % determine if we're doing a phase space calculation or not
    if(~isempty(dataPDE))
        phasespace = dataPDE.opts.phasespace;
    else
        phasespace = dataSDE.opts.phasespace;
    end

    % determine which results we have and extract the relevant information 
    if(~isempty(dataPDE))
        doPDE = true;
        rho = dataPDE.rho_t;
        plotTimes = dataPDE.plotTimes;

        if(~phasespace)
            qLine = dataPDE.aLine;
        else
            qLine = dataPDE.qLine;
            pLine = dataPDE.pLine;
            phi = dataPDE.phi_t;
        end
    else
        doPDE = false;
    end

    if(~isempty(dataSDE))
        doSDE = true;
        if(~doPDE)
            plotTimes = dataSDE.tStoc;
        end
    else
        doSDE = false;
    end

    if(nargin<3)
        optsPlot = struct;
    end
    optsPlot.doPlots = true; % force plotting on
    
    % set up a figure
    if(isfield(optsPlot,'hf'))
        hf = optsPlot.hf;
        haq = optsPlot.haq;
        if(phasespace)
            hap = optsPlot.hap;
        end
    else
        if(phasespace)
            hf = figure('Position',[100,100,1200,800]);
            haq = subplot(1,2,1);
            hap = subplot(1,2,2);
            optsPlot.hap = hap;
        else
            hf = figure('Position',[100,100,800,800]);
            haq = axes;
        end
        optsPlot.haq = haq;
    end
    
    nPlots = length(plotTimes);
    
    optsPlot.xlabel = '$\theta$';
    optsPlot.ylabel = 'Density';
    
    optsPlot.plain = true;
    optsPlot.linewidth = 1.5;
    linecolors=linspecer(nPlots);

    % set up the plotting 
    optsPlotq = optsPlot;
    optsPlotq.ha = optsPlot.haq;
    if(phasespace)
        optsPlotp = optsPlot;
        optsPlotp.ha = optsPlot.hap;
    end
    
    % whether to shift the CoM of the SDE calculation (really only useful
    % for no external potential)
    
    if(doSDE && isfield(optsPlot,'shiftCoM'))
        shiftCoM = optsPlot.shiftCoM;
    elseif(doSDE)
        if(strcmp(dataSDE.opts.optsVext.VextFn,'VextZero'))
            shiftCoM = true;
        else
            shiftCoM = false;
        end
    else
        shiftCoM = false;
    end
    
    if(shiftCoM)
        qLine = dataPDE.aLine;
        rhoFinal = dataPDE.rho_t(:,end);
        theta = dataPDE.aLine.Pts.y;
        Int = dataPDE.aLine.Int;
        
        opts.targetCoM = Int*(theta.*rhoFinal);
        xShifted = KuramotoShiftSDE(dataSDE.x,opts);
        dataSDE.xRaw = dataSDE.x;
        dataSDE.x = xShifted;
    end
    
    dataSDEq = dataSDE;
    if(phasespace)
        dataSDEp = dataSDE;
        dataSDEp.x = dataSDEp.p;
        if(doPDE)
            optsPlotp.xLim = [dataPDE.pLine.yMin,dataPDE.pLine.yMax];
        else
            minp = min(dataSDEp.p(:));
            maxp = max(dataSDEp.p(:));
            optsPlotp.xLim = [minp,maxp];
        end
        optsPlotp.xLabel = '$p$';
    end
    
    % plot the results
    
    setYLimRho = true;
    if(isfield(optsPlot,'rhoMax'))
        rhoMax = optsPlot.rhoMax;
    elseif(doPDE)
        rhoMax = 1.1*max(rho(:));
    else
        setYLimRho = false;
    end
    
    setYLimPhi = true;
    if(phasespace)
        if(isfield(optsPlot,'phiMax'))
            phiMax = optsPlot.phiMax;
        elseif(doPDE)
            phiMax = 1.1*max(phi(:));
        else
            setYLimPhi = false;
        end 
    end

    
    for iPlot = 1:nPlots
        hold(haq,'on');
        
        linecolor = linecolors(iPlot,:);
        displayname = ['$t = ' num2str(plotTimes(iPlot)) '$'];
        
        if(doPDE)
            optsPlotq.handlevisibility = 'on';
            optsPlotq.linecolor = linecolor;
            optsPlotq.displayname = displayname;
            qLine.plot(rho(:,iPlot),optsPlotq);
           
            if(phasespace)
                hold(hap,'on');
                optsPlotp.handlevisibility = 'on';
                optsPlotp.linecolor = linecolor;
                pLine.plot(phi(:,iPlot),optsPlotp);
            end
        end
                        
        if(doSDE)
            if(doPDE)
               optsPlotq.handlevisibility = 'off';
            end
            optsPlotq.tPlot = plotTimes(iPlot);
            optsPlotq.linecolor = linecolor;
            PlotKuramotoSDE(dataSDEq,optsPlotq);

            if(phasespace)
                if(doPDE)
                    optsPlotp.handlevisibility = 'off';
                end
                optsPlotp.tPlot = plotTimes(iPlot);
                optsPlotp.linecolor = linecolor;
                PlotKuramotoSDE(dataSDEp,optsPlotp);
            end
        end
        
        if(setYLimRho)
            ylim(haq,[0,rhoMax]);
        end
        
        if(phasespace && setYLimPhi)
            ylim(hap,[0,phiMax]);
        end

    end   
    
    legend('interpreter','latex');
    
    if(isfield(optsPlot,'doPDF') && optsPlot.doPDF)
        save2pdf(optsPlot.snapshotsFile,hf);
    end

end