%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Green's function simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  clear all, clc, 
%  close all
%-------------------------------------------------
% Constants and parameters
%-------------------------------------------------
c = 343;
lx = 1;
ly = 1;

%source position
Ps = [ 0.5, 0.5 ];
Pr = [ 0.2 0.2 0.2];    %receiver position

% Green function frequencies
f = 20:1:1200;       % resolution of 10Hz
k = 2*pi*f/c;       % frequency dependent constant k
alpha = 0.2;
V = lx*ly;
T_60=0.16*V/(alpha*(lx*ly)); % sabine
tau = T_60/13.8;

% ----------------------------------------------------------------------
% Calculate Greens function
% ----------------------------------------------------------------------
% number of modes in each direction
Nx = ceil(max(k)*lx/pi);
Ny = ceil(max(k)*ly/pi);

fprintf('\nNumber of modes: %i, %i\n', Nx, Ny)
fprintf('Frequency: %f\n\n Hz', max(f))

P = greens.greensNeumannFreq2D(Nx,Ny,lx,ly,Ps,Pr,k,c,tau);

%-------------------------------------------------
% PLOT simulated Greens function
%-------------------------------------------------
figure()

G_sim_dB = 20*log10(abs(P));
semilogx(f,G_sim_dB,'r');



