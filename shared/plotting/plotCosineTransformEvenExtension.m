function plotCosineTransformEvenExtension(x, boundaryCond, fs, ax, titlefig) % 'Gaussian cosine transform' ax3 = nexttile;    
    X = dctmodes(x, boundaryCond);
    NDCT = length(X); % DCT resolution
    freqs = fs/2*(0:NDCT-1)/(NDCT);
    
    stem(ax, freqs, X)
    title(ax,titlefig)
    xlabel(ax,'Frequency (Hz)')   
    ylabel(ax,'Magnitude');
    grid
end