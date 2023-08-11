function [domain] = domainTestSetup(nmodes)
    xminmax = [0,10];
    l = xminmax(2) - xminmax(1);
    c = 343;
    rho = 1.2;
    
    wavelength_min = l/nmodes;
    fn = c/(2*wavelength_min); % Nyquist frequency
    fs = 2*fn;                 % Sampling frequency
    
    dx = wavelength_min/2;
    dt = 1/fs;
    
    x1d = xminmax(1):dx:xminmax(2);
    
    domain = Domain1D(x1d,dx,dt,fs,c,rho,xminmax,BoundaryType.Neumann);
end