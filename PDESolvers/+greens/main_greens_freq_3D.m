%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Green's function simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  clear all, clc, 
%  close all
 
%-------------------------------------------------
% Constants and parameters
%-------------------------------------------------
c = 343;
lx = 4;  
ly = 2.7;      
lz = 3;
V  = lx*ly*lz;

Ps_all=[3, 1.5, 2 ];                %source position 
Pr_all=[ 0.5 1.7 0.5;1.0 1.7 1];    %receiver position
% Ps_all=[0.5, 0.5, 0.5 ];                %source position 
% Pr_all=[ 0.1 0.1 0.1];    %receiver position

alpha_all = [0.3227];

% reverberation time 
for aind=1:length(alpha_all)
    
    alpha=alpha_all(aind);
    T_60=0.16*V/alpha/2/(lx*ly+lz*ly+lx*lz);
    tau = T_60/13.8;

    % Green function frequencies
    
    f = 20:1:1000;       % resolution of 10Hz
    k = 2*pi*f/c;       % frequency dependent constant k

    % ----------------------------------------------------------------------
    % Calculate Greens function
    % ----------------------------------------------------------------------
    % number of modes in each direction
    Nx = ceil(max(k)*lx/pi);
    Ny = ceil(max(k)*ly/pi);
    Nz = ceil(max(k)*lz/pi);
    
    fprintf('\nNumber of modes: %i, %i, %i\n', Nx, Ny, Nz)
    fprintf('Frequency: %f\n\n Hz', max(f))

    P = greens.greensNeumannFreq3D(Nx,Ny,Nz,lx,ly,lz,Ps_all,Pr_all,tau,k,c);
    %[P, f] = greens.greensDirichletFreq3D(Nx,Ny,Nz,lx,ly,lz,Ps_all,Pr_all,c);

    %-------------------------------------------------
    % PLOT simulated Greens function
    %-------------------------------------------------
    %%
    figure()
    G_sim_dB = 20*log10(abs(P));
    %semilogx(f,G_sim_dB,'r');
    semilogx(f,G_sim_dB,'r');
end



