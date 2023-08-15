clear all
close all

import setup.*
import models.*
import models.types.*

do_animation = false;
write_gif = do_animation && false;

write_data = true;
write_fig = false;

base_path = '~/data/fm-sem/';

%% SETUP SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%tmax = 0.0135;
tmax = 0.2;

% .FDTD_FDTD | .FOURIER_FOURIER | .FOURIER_FDTD | .SEM_SEM | .FOURIER_SEM
coupling_methods = "FOURIER_SEM"; 
source_partition = "LEFT";
x0_pos = 2.5;
frac_partition1 = 0.50;

scheme_order = 6; % 2 or 6

% SEM SETUP
P_order = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xminmax = [0,10]; % domain length
fmax = 1000;      % [Hz]
c = 343;          % m/sec speed
rho = 1.2;

interface_order = scheme_order;

%% SETUP and SOLVE
switch coupling_methods
    case "FDTD_FDTD"
        ppw1 = 4;
        ppw2 = 4;

        [dx1] = utilsDD.calcSpatial(xminmax(2),fmax,ppw1,c);
        [dx2] = utilsDD.calcSpatial(xminmax(2),fmax,ppw2,c);
        
        [p1,p2,s1,s2,iter] = run_FDTD_FDTD(tmax, xminmax, x0_pos, dx1, dx2, fmax, c, rho, scheme_order);
    case "FOURIER_FOURIER"
        % ppw are restricted by the interface scheme
        ppw1 = 4;
        ppw2 = 12;
        
        [dx1] = utilsDD.calcSpatial(xminmax(2),fmax,ppw1,c);
        [dx2] = utilsDD.calcSpatial(xminmax(2),fmax,ppw2,c);
        
        [p1,p2,s1,s2,iter] = run_FOURIER_FOURIER(tmax, xminmax, x0_pos, dx1, dx2, fmax, ...
            source_partition, c, rho, scheme_order);
    case "FOURIER_FDTD"
        ppw1 = 4;
        ppw2 = 8;
        
        [dx1] = utilsDD.calcSpatial(xminmax(2),fmax,ppw1,c);
        [dx2] = utilsDD.calcSpatial(xminmax(2),fmax,ppw2,c);
        
        [p1,p2,s1,s2,iter] = run_FOURIER_FDTD(tmax, xminmax, x0_pos, dx1, dx2, fmax, c, rho, scheme_order);
    case "SEM_SEM"
        % TODO: assert failing: equal number of elements required

        ppw1 = 2;
        ppw2 = 4;
        
        [dx1] = utilsDD.calcSpatial(xminmax(2),fmax,ppw1,c);
        [dx2] = utilsDD.calcSpatial(xminmax(2),fmax,ppw2,c);
        dx_i = dx2;
        iface = SEMUniformInterface(dx_i, 3, 1, false, InterfaceLocation1D.LEFT);
        
        [p1,p2,s1,s2,iter] = run_SEM_SEM(tmax, xminmax, x0_pos, dx1, dx2, ...
            fmax, c, rho, P_order, interface_order, iface);    
    case "FOURIER_SEM"                
        ppw1 = 4;
        ppw2 = 4;
        
        xminmax_1 = xminmax;
        xminmax_2 = xminmax;
        
        xminmax_1(2) = round(xminmax(2)*frac_partition1);
        xminmax_2(2) = xminmax(2) - xminmax_1(2);
        
        [dx1] = utilsDD.calcSpatial(xminmax_1(2),fmax,ppw1,c);
        dx2 = dx1*ppw1/ppw2;
        %[dx2] = utilsDD.calcSpatial(xminmax_2(2),fmax,ppw2,c);
        dx_i = dx2/2;
        iface = SEMUniformInterface(dx_i, 3, 1, false, InterfaceLocation1D.LEFT);
        
        % adjust SEM length for integer elements
        xminmax_2(2) = ceil(xminmax_2(2)/(dx2*P_order))*(dx2*P_order);
        fprintf('SEM domain length: %0.4f\n', xminmax_2(2))
        
        if iface.Nk > 0
            if iface.loc == InterfaceLocation1D.LEFT
                xminmax_2(1) = xminmax_2(1)-iface.dx*iface.Nk;
            else
                xminmax_2(2) = xminmax_2(2)+iface.dx*iface.Nk;
            end
        end
        
        [p1,p2,s1,s2,iter] = run_FOURIER_SEM(tmax, xminmax_1, xminmax_2, ...\
            dx1, dx2, fmax, source_partition, x0_pos, c, rho, P_order, interface_order, iface);
    otherwise
        error('coupling method not supported');
end

x1d_1 = s1.domain.x1d;
x1d_2 = s1.domain.x1d(end) - s2.domain.x1d(1) + s2.domain.x1d;

% normalize
tnorm = 0.007;
n_norm = round(tnorm/s1.domain.dt);
if source_partition == "LEFT"
    norm_factor = max(p1(n_norm,:));
elseif source_partition == "RIGHT"
    norm_factor = max(p2(n_norm,:));
else
    error('no such source location')
end

p1 = p1/norm_factor;
p2 = p2/norm_factor;

% ref
dx1 = s1.domain.dx; dx2 = s2.domain.dx;
dt1 = s1.domain.dt; dt2 = s2.domain.dt;
dt = max(dt1,dt2);

if source_partition == "LEFT"        
    x0_adj = s1.src.x0;
elseif source_partition == "RIGHT"
    x0_adj = s2.src.x0;
else
    error('no such source location')
end

xminmax_full = [x1d_1(1),x1d_2(end)];

[p_ref, s_ref] = calcReference(xminmax_full, ppw1, dt, fmax, x0_adj, c, rho, iter, tnorm);

if write_data
    L = xminmax_full(2) - xminmax_full(1);
    save(sprintf('%s/%s_ppw%i_%i_src_%s_fmax_%i.mat', ...\
        base_path, coupling_methods, ppw1, ppw2, source_partition, fmax), ...\
        'p1', 'p2', 'p_ref', 's1', 's2', 's_ref', 'ppw1', 'ppw2', 'L', 'interface_order', 'fmax', 'iter');
end

%% PLOT
title_str = sprintf('%s-%s \n dx ratio 1 : %i   dt ratio %i : 1', .../
    s1.solver_type, s2.solver_type, ppw2/ppw1, s1.domain.dt/s2.domain.dt);
legend1 = sprintf('%s',s1.solver_type);
legend2 = sprintf('%s',s2.solver_type);

if write_fig
    h = figure(1);
    set(h,'Position',[10 10 1600 800])
    tiledlayout(2,1)
    ax1 = nexttile;
    ax2 = nexttile;
    
    i1 = utilsDD.getClosestIndex(x1d_1, 2.5);
    i2 = utilsDD.getClosestIndex(x1d_1, 5.0);

    err_dB_iface = max(20*log10(abs(p1(iter,i1:i2))));

    fprintf('%s_ppw%i_%i_fmax_%i\n',coupling_methods,ppw1,ppw2,fmax);
    fprintf('Err interface: %i dB\n', round(err_dB_iface))

    plotting.plotSemiLogCoupled(x1d_1, x1d_2, s_ref.domain.x1d, p1(iter,:), p2(iter,:), p_ref(iter,:), title_str, legend1, legend2, ax1, ax2)

    figure(3); 
    plot(x1d_1(i1:i2), 20*log10(abs(p1(iter,i1:i2))))
end

if do_animation
    h = figure(2);
    set(h,'Position',[10 10 1600 800])
    tiledlayout(2,1)
    ax1 = nexttile;
    ax2 = nexttile;
    
    for n=1:iter
        if mod(n-1,5) == 0
            t = n*s1.domain.dt;
            fprintf('n=%i, t=%0.4f\n',n,t)
            plotting.plotSemiLogCoupled(x1d_1, x1d_2, s_ref.domain.x1d, p1(n,:), p2(n,:), p_ref(n,:), ...\
                title_str, legend1, legend2, ax1, ax2)
            drawnow

            if write_gif
                path = sprintf('%s/%s_ppw%i_%i_fmax_%i',base_path,coupling_methods,ppw1,ppw2,fmax);
                plotting.writeGif(h, path, n==1);
            end
        end
    end    
end

function [p, s1] = calcReference(xminmax, ppw, dt, fmax, x0_pos, c, rho, iter, tnorm)
    import models.types.*

    [dx,nmodes] = utilsDD.calcSpatial(xminmax(2),fmax,ppw,c);
    if dt > dx/c
        error('dt not satisfying CFL!')
    end
    
    x1d = xminmax(1):dx:xminmax(end);
    assert(abs(x1d(end) - xminmax(end)) < eps('double'))
    
    sourceFactor = utilsDD.sourceScalingFourier()*2;
    
    domain1 = models.Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = utilsDD.defaultSource(fmax,sourceFactor,domain1,x0_pos);
    assert(abs(src1.x0 - x0_pos) < eps('single'))
    
    s1 = models.Simulation1D(SolverType.FOURIER, domain1, src1, solvers.CustomFourier(nmodes));
    
    p = solvers.runSingleDomainSolver(iter, s1);
    
    n_norm = round(tnorm/dt);
    norm_factor = max(p(n_norm,:));
    
    p = p/norm_factor;
end