function plotTimeSignal(ts,x,ax,titleFig)
    plot(ax,ts,x)
    title(ax,titleFig)
    xlabel(ax,'time [s]')    
    ylabel(ax,'Amplitude');
    grid
end