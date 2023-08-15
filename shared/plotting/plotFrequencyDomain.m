function [fROI, XXROI] = plotFrequencyDomain(ax, x, fs, f_pulse, from, title_str, legend_style, spectrum_type)    
    if nargin == 7
        spectrum_type = "SPL";        
    end
    
    NFFT = length(x); %8*1024; % FFT resolution
    
    [fROI, XXROI] = transformToFreq(x,fs,f_pulse,from,spectrum_type,NFFT);
    
    if spectrum_type == "POWER"
        ylabel_str = 'Power Spectrum [dB]';        
        semilogx(ax,fROI,10*log10(XXROI),legend_style, 'LineWidth', 1.5);        
    elseif spectrum_type == "SPL"
        ylabel_str = 'SPL [dB]';
        semilogx(ax,fROI,20*log10(abs(XXROI)),legend_style, 'LineWidth', 1.5);
    else
        error('Spectrum type not supported')
    end
    
    title(ax,title_str);
    %xlabels = [from, Fs_pulse/8, Fs_pulse/4, Fs_pulse/2];
    %set(gca,'Xtick',xlabels)
    %set(gca,'XtickLabel',xlabels)
    xlabel(ax,'Frequency [Hz]')    
    ylabel(ax,ylabel_str, 'Interpreter','latex');
    xlim(ax,[from,f_pulse])
    grid(ax,'on')
end
        
function [fROI, XXROI] = transformToFreq(x, fs, f_upper, from, spectrum_type,  NFFT)
    N = length(x);    

    freqs = fs*(0:NFFT/2-1)/NFFT;

    if spectrum_type == "POWER"
        X = fftshift(fft(x, NFFT));
        XX = X.*conj(X)/(N*N); %computing power with proper scaling        
    elseif spectrum_type == "SPL"
        XX = fftshift(fft(x, NFFT)/NFFT);        
    else
        error('Spectrum type not supported')
    end

    XXHalf = XX(NFFT/2+1:end); % multiply with 2?

    indxFrom = 1+floor(from/freqs(end) * length(freqs));
    indxTo = min(floor(f_upper/freqs(end) * length(freqs)),NFFT/2);
    fROI = freqs(indxFrom:indxTo);
    XXROI = XXHalf(indxFrom:indxTo);
end