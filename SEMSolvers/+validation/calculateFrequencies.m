function [Pfreqs, XY, freq_bins] = calculateFrequencies(solver_f,c,rho,freq_min_max,t)    
    % relation between frequency and time
    fs = 2*freq_min_max(2);
    N = ceil(fs*t);
    fbin_res = fs/N;
    
    freq_bins = flip(freq_min_max(2):-fbin_res:freq_min_max(1));
    
    L = length(freq_bins);
    
    i = 1;
    
    % first iteration outside loop to determine number of grid points
    params = models.AcousticParameters(freq_bins(i),c,rho);        
    fprintf('Solving for frequency %0.2f...\n', params.f)        
    [Pfreq, XY] = solver_f(params.k);
    
    Pfreqs = zeros(L,size(XY,1));
    Pfreqs(i,:) = Pfreq(:);    
    
    XY_prev = XY; 
    
    for i=2:L
        params = models.AcousticParameters(freq_bins(i),c,rho);
        
        fprintf('Solving for frequency %0.2f...\n', params.f)        
        [Pfreq, XY] = solver_f(params.k);
                
        assert(all(XY_prev == XY, 'all')); XY_prev = XY; % XY should be the same for all simulations
        
        Pfreqs(i,:) = Pfreq(:);
    end
end