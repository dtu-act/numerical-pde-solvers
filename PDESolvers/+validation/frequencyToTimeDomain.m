function [X_time] = frequencyToTimeDomain(X_freqs,padding)
%
% INPUT: X_freqs is a M x N matrix with
%           row   : frequency index
%           column: node index
%
% OUTPUT: time_solutions being a T X N matrix with
%           row   : time index
%           column: node index
%
% EXAMPLE: X_time(i,j) refers to the ith frequency amplitude at node j
%
%
    M = size(X_freqs,1); % number of frequencies
    
    if nargin == 1
        T = (M-1)*2;         % number of time steps        
    else
        T = padding;         % number of time steps
    end    
            
    N = size(X_freqs,2); % number of nodes
    
    X_time = zeros(T,N);    
    tf = zeros(M,1);
    
    for j=1:N
        tf(:) = X_freqs(:,j);
        
        % we need double-sided: make it even and symmetric, i.e. time
        % domain signal will be pure real        
        padding = T - (length(tf) + length(tf) - 2);
        tf_double = [tf; zeros(padding,1); flip(conj(tf(2:end-1)))];
        
        ir = ifft(tf_double); % real(ifft(tf_double, padding))
        
        assert(all(imag(ir) == 0, 'all'))
        assert(T == length(ir))
        
        X_time(:,j) = ir(:);
    end    
end

% https://se.mathworks.com/help/matlab/ref/fft.html
% https://se.mathworks.com/matlabcentral/answers/109444-fft-ifft-single-and-double-spectrum
% https://electronics.stackexchange.com/questions/12407/what-is-the-relation-between-fft-length-and-frequency-resolution
% https://se.mathworks.com/matlabcentral/answers/278387-fft-conversion-of-sound-pressure-in-the-time-domain
