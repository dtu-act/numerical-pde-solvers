clear all
close all
addpath('../shared')

% SIMULATION PARAMETERS
c = 343.0;      % m/s speed of sound
Porder = 4;

%[20 40 80 160];
%M = 500;         % elements
L = 2.5;
source_type = "gaussian"; % "point", "mmf", "gaussian"

f = 300; % Hz
fmax = 1000;
ppw = 10;

k = f*2*pi/c;
N_wavefun_greens = 5000;
x0 = L/2;



switch source_type
    case 'mmf'
        ufun = @(x) sin(pi*x);
        F = @(x) -pi^2*sin(pi*x) + k^2*sin(pi*x);
        source = models.SourceModel(models.SourceType.Function, F);
        bounds = [ufun(0), ufun(L)];
    case 'point'
        source = models.PointSource(x0,-1);
        bounds = [0,0];
    case 'gaussian'
        sigma = c/(pi*fmax/2);
        F_gauss = @(x) -(exp(-(x-x0).^2/(sigma^2)));
        source = models.SourceModel(models.SourceType.Function,F_gauss);
        bounds = [0,0];
    otherwise
       error('source_type not supported')
end

[ps, x1d] = solvers.solveHelmholtz1D(L, k, Porder, c, fmax, ppw, source, bounds);

switch source_type
    case "mmf"
        plot(x1d, ps)        
        err_l2 = norm( ps - ufun(x1d) );
        fprintf("Err l2: %f\n", err_l2)
        hold on
        plot(x1d, ufun(x1d), '--o')
        legend('SEM','Analytical')
        hold off
    case "point"        
        Np = length(ps);

        disp("Calc greens function....")
        Pref  = greens.greensDirichletSpatial1D(N_wavefun_greens,x1d',x0,k);

        plot(x1d, ps, '-o')

        hold on
        plot(x1d, Pref, '--.')
        hold off
        legend('SEM 1D','Greens 1D')

        err_l2 = norm( ps - Pref, 2); %/Np;
        fprintf("Err l2: %f\n", err_l2)
    case "gaussian"
        Np = length(ps);

        disp("Calc greens function....")
        Pref = greens.greensDirichletSpatial1D(N_wavefun_greens,x1d,x0,k);

        r = 0.4;
        fd = @(x) abs(x-x0) -  r;
        psim_roi = validation.extractROI(ps,x1d,fd);
        Pref_roi = validation.extractROI(Pref,x1d,fd);
        C = validation.calcNormConstant(psim_roi,Pref_roi);
        
        plot(x1d, C*psim_roi, '-o')

        hold on
        plot(x1d, Pref_roi, '--o')
        hold off
        legend('SEM 1D','Greens 1D')

        err_l2 = norm( C*psim_roi - Pref_roi, 2); %/Np;
        fprintf("Err l2: %f\n", err_l2)
    otherwise
        error('source_type not supported')
end