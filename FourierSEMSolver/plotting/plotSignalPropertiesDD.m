function [] = plotSignalPropertiesDD(ts, x, boundaryCond, fs, fs_pulse, fromHz, dct_type, fig_i)
    figure(fig_i)
    tiledlayout(4,1) % Requires R2019b or later

    ax1 = nexttile;
    plotTimeSignal(ts,x,ax1,'Time domain');

    ax2 = nexttile;
    plotPowerSingleSided(x,fs,fs_pulse,fromHz,ax2,'Power Spectral Density');

    ax3 = nexttile;
    plotCosineTransform(x,fs,dct_type,ax3,'Cosine transform')
    
    ax4 = nexttile;
    plotCosineTransformEvenExtension(x,boundaryCond,fs,ax4,'FFT transform with even extension')
end