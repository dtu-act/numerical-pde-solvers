function [vx, etov, Ms] = setupUniformGrid1D(c,fmax,P,ppw,xl,xu)
    wavelength_min = c/fmax;
    dx = wavelength_min/ppw;
    Ns = length(linspace(xl,xu,ceil((xu-xl)/dx)));

    Ms = ceil(Ns/P);
    vx = linspace(xl,xu,Ms+1);
    
    etov = sem.setupEToV1D(vx);
    %vxconn = 1:length(etov);
end