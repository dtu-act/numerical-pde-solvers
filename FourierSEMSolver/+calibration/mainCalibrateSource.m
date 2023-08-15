clear all
close all

%% SETUP SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

method = SolverType.SEM;
order = 4;
iter = 250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_source = 1000;  % [Hz]
nmodes = 100;       % number of modes to include
xminmax = [0,5];    % domain length    
c = 343;            % m/sec speed
rho = 1.2;

%% MESH
L = xminmax(2)-xminmax(1);
wavelength_min = L/nmodes;
fn = c/(2*wavelength_min); % Nyquist frequency
fs = 2*fn;                 % Sampling frequency

% spatial/temporal relation: dx = c/(fmax*2); dt = 1/(fmax*2);
%dx = wavelength_min/2; % == c/Fs;

x1d = linspace(xminmax(1),xminmax(2),nmodes*2 + 1);
dx = x1d(2)-x1d(1);

sourceFactor = calibrateSource(method, dx, nmodes, c, rho, xminmax, fn_source, order, iter);

fprintf('Source factor: %.20fe\n', sourceFactor);

