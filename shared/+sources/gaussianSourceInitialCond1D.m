function [source, source_x, mu] = gaussianSourceInitialCond1D(c, fmax)
    sigma_t = 1/(pi*fmax/2);
    sigma_x = c/(pi*fmax/2);
    mu = 0; % no smooth initialization when used as initial condition
    
    source_t = @(t) exp(-(t-mu).^2/(sigma_t^2));        
    source_x = @(x,t,x0) 0.5*exp(-(x-x0 - c*t).^2/(sigma_x^2)) + 0.5*exp(-(x-x0 + c*t).^2/(sigma_x^2));
    
    source = @(x,t,x0) source_t(t).*source_x(x,0,x0);
end