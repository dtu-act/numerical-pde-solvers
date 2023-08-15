function plotCosineTransform(x, fs, dct_type, ax, titlefig) % 'Gaussian cosine transform' ax3 = nexttile;
    NDCT = length(x); % DCT resolution
    freqs = fs/2*(0:NDCT/2-1)/(NDCT/2);
    X = dct(x,NDCT,'Type',dct_type); 
    stem(ax, freqs, X(1:NDCT/2))
    title(ax,titlefig)
    xlabel(ax,'Frequency (Hz)')   
    ylabel(ax,'Magnitude');
    grid
end