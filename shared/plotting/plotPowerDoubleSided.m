function plotPowerDoubleSided(x, L, Fs, ax, title)
    NFFT = 1024; % FFT resolution
    %Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
    X = fftshift(fft(x, NFFT));
    
    freqs = Fs*(-NFFT/2:NFFT/2 - 1)/NFFT; %linspace(0,Fs/2,L);
    Pxx=X.*conj(X)/(L*L); %computing power with proper scaling
    %semilogx(freqs,10*log10(Pxx),'r');
    plot(freqs,10*log10(Pxx),'r');
    %set(gca, 'XScale', 'log')
    %title(ax,title);
    set(gca,'XTickLabel',arrayfun(@(x) num2str(x), get(gca,'XTick'),'un',0))
    xlabel(ax,'Frequency (Hz)')    
    ylabel(ax,'Power Spectral Density- P_{xx} dB/Hz');
    grid
end

% https://se.mathworks.com/help/matlab/ref/matlab.graphics.axis.decorator.numericruler-properties.html