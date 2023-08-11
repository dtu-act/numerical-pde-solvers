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

% Green function frequencies

from_fundamental_freq = false;

if from_fundamental_freq    
    % set number of wave functions and calculate frequency from this
    Nx = 4;
    
    k = pi*Nx/lx;
    Ny = ceil(k*ly/pi);
    
    lambda = 2*pi/k;
    f = c/lambda;
else
    % set frequency and calculate number of wave functions needed from this
    f = 500;
    
    lambda = c/f;    
    k = 2*pi/lambda;
    
    % number of modes in each direction
    Nx = 50; %ceil(k*lx/pi);
    Ny = 50; %ceil(k*ly/pi);
end

fprintf('\nNumber of modes: %i, %i\n', Nx, Ny)
fprintf('Frequency: %f Hz\n\n', f)

% ----------------------------------------------------------------------
% Calculate Greens function
% ----------------------------------------------------------------------

%[P, xs, ys] = greens.greensNeumannSpatial2D(Nx,Ny,lx,ly,Ps,k);
[P, X, Y] = greens.greensDirichletSpatial2D(Nx,Ny,lx,ly,Ps,k);

%-------------------------------------------------
% PLOT simulated Greens function
%-------------------------------------------------
plotting.plotHelmholtz(abs(P), X, Y, 'Greens')

% figure(5)
% contourf(xs, ys, P);
% colorbar;