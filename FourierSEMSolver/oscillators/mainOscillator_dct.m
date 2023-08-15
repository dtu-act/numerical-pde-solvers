%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

import models.*
import models.types.*

boundarycond = BoundaryType.Neumann;
%boundarycond = BoundaryType.Dirichlet;

%load_signal = "none";
%load_signal = "test-data/pressures_even_source";
%load_signal = "test-data/pressures_odd_source";
%load_signal = "test-data/pressures_even_boundary";
load_signal = "test-data/pressures_even_neumann";
%load_signal = "test-data/pressures_odd_boundary";

if load_signal ~= "none"
    data = load(load_signal, 'p_current', 'x1d', 'fn', 'fs', 'length');
    x = data.p_current;    
    x1d = data.x1d;        
    l = data.length;
    fn = data.fn;
    fs = data.fs;
else   
    x1d = linspace(0,2*pi,100);
    x = sin(x1d);
    fn = 1;
    l = 5;
end

if boundarycond == BoundaryType.Dirichlet
    % 1 - convert a time signal to frequency domain
    X = dct.dctmodes([0,x,0], boundarycond);
    % 2 - convert back from frequency domain to time domain and compare with original
    x_idct = dct.idctmodes(X, boundarycond);
    x_idct = x_idct(2:end-1);
else
    % 1 - convert a time signal to frequency domain
    X = dct.dctmodes(x, boundarycond);    
    % 2 - convert back from frequency domain to time domain and compare with original
    x_idct = dct.idctmodes(X, boundarycond);
end

N = length(x);
NX = length(X);
freqNyquest = (0:NX-1)/(NX-1)*fn;

figure(2)
stem(freqNyquest,X)

switch boundarycond
    case BoundaryType.Neumann
        Xdouble = [X,fliplr(X(2:end-1))];
        x_osc = idctL(x1d,Xdouble,l);
        
        figure(3)        
        a1 = plot(x1d,x, 'x--r');
        hold on
        a2 = plot(x1d,x_osc(1:N), '-g');        
        a3 = plot(x1d,x_idct, '-b');
        %a4 = plot(x1d,MinvN(1:N), '*-m');
        xlabel('length [m]')
        ylabel('Magnitude');
        legend([a1,a2,a3], ["Original signal", "Oscillator transformed", "iDCT"]);
        hold off
    case BoundaryType.Dirichlet
        figure(3)        
        a1 = plot(x1d,x, 'x--r');
        hold on
        a2 = plot(x1d,x_idct, '-b');
        xlabel('length [m]')
        ylabel('Magnitude');
        legend([a1,a2], ["Original signal", "iDCT"]);
        hold off
    otherwise
        error('Boundary condition not supported')
end

function [p] = idctL(x,m,l)
    N = length(m);
    
    p = 0;
    for n=1:N
        p = p + m(n)*cos(pi/l*(n-1)*x);
    end
end

function [p] = idctN(m)
    N = length(m);

    ks = 0:(N-1);
    p = 0;

    for n=1:N
        p = p + m(n)*cos(2*pi*ks*(n-1)/N);
    end
end
