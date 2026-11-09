function KuramotoMovie(dataPDE,dataSDE,optsPlot)

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
    
    % determine if we're writing an mp4
    doMp4 = (isfield(optsPlot,'doMp4') && optsPlot.doMp4);
    if(doMp4)
        if(isfield(optsPlot,'mp4File'))
            v = VideoWriter(optsPlot.mp4File,'MPEG-4');
        else
            v = VideoWriter('movieDemo.mp4','MPEG-4');
            optsPlot.mp4File = 'movieDemo.mp4';
        end
        v.FrameRate=10;
        v.Quality = 100;
        open(v)
    end

    
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
    if(isfield(optsPlot,'linecolor'))
        optsPlot.linecolor = optsPlot.linecolor;
    else
        optsPlot.linecolor = 'r';
    end
    
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
        hold(haq,'off');
        
        if(doPDE)
            qLine.plot(rho(:,iPlot),optsPlotq);
           
            if(phasespace)
                hold(hap,'off');
                pLine.plot(phi(:,iPlot),optsPlotp);
            end
        end
                
        if(doPDE && doSDE)
            hold(haq,'on');
            if(phasespace)
                hold(hap,'on');
            end
        end
        
        if(doSDE)
            optsPlotq.tPlot = plotTimes(iPlot);
            PlotKuramotoSDE(dataSDEq,optsPlotq);

            if(phasespace)
                optsPlotp.tPlot = plotTimes(iPlot);
                PlotKuramotoSDE(dataSDEp,optsPlotp);
            end
        end
        
        if(setYLimRho)
            ylim(haq,[0,rhoMax]);
        end
        
        if(phasespace && setYLimPhi)
            ylim(hap,[0,phiMax]);
        end

        if(doMp4)
            frame = getframe(hf);
            writeVideo(v,frame);
        end
        
        pause(0.2)
%         pause
    end

    if(doMp4)
        fprintf(['Movie saved in ' optsPlot.mp4File '\n']);
        close(v);
    end
    
    

end