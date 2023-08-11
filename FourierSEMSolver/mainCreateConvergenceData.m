clear all
close all

write_data = false;
write_gif = true;
do_plot = true;
base_path = '~/data/fm-sem/dxFourier_x2_dxInterface_x1';

%% SETUP SIMULATION PARAMETERS

% "fourier" | "sem" | "fourier_sem"
conv_test = "fourier_sem";

refine_factor = 4; % 1 = Nyquist

ppw_sem = 2*refine_factor;  % 2 corresponds to Nyquist
ppw_fourier = max(ppw_sem/2,2);

scheme_order = 4;
interface_order = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmax = 0.2; %0.0115/1.4; % [sec]
x0_pos = 1.25;

fmax = 1000;  % Nyquist frequency [Hz]
c=343;
rho=1.2;

% spatial resolutions
l_d1 = 5/2;
l_d2 = l_d1;

dx1 = calcSpatial(l_d1,fmax,ppw_fourier,c);
dx2_i = calcSpatial(l_d1,fmax,ppw_sem,c);
iface = SEMUniformInterface(dx2_i, 3, 1, false, InterfaceLocation1D.LEFT);

l_d2_inner = l_d1 - iface.dx*iface.Nk;

[dx2,~] = calcResolutionEvenElements(l_d2_inner,fmax,ppw_sem,scheme_order,c);

cfl = sqrt(3);
    
if conv_test == "fourier"
    dt = dx1/(c*cfl);
    iter = tmax/dt;
    if rem(iter,1) ~=0
        iter = ceil(iter);
        dt_adjust = tmax/iter;
        fprintf('WARNING: dt1 adjusted from %e to %e to match tmax', dt, dt_adjust)
        dt = dt_adjust;
    end
    
    L = l_d1+l_d2;
    xminmax = [0,L];
    
    % round to match number of modes (discretization will differ from SEM)
    %[dx,nmodes] = calcSpatial(L,fmax,ppw_fourier,c);
    nmodes = L/(2*dx1);
    x1d = xminmax(1):dx1:xminmax(end);
    
    sourceFactor = sourceScalingFourier();
    
    domain1 = Domain1D(x1d, dx1, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fmax,sourceFactor,domain1,x0_pos);
    s1 = Simulation1D(SolverType.FOURIER, domain1, src1, CustomFourier(nmodes));
    
    p = runSingleDomainSolver(iter, s1);
    p = p*2;
    
    dx=dx1;
    
    if write_data
        save(sprintf('%s/%s_L%0.2f_dx%0.3f_dt%0.6f.mat', base_path, conv_test, L, dx1, dt), ...\
            'p', 'x1d', 'dx', 'dt', 'tmax', 'L');
    end
elseif conv_test == "sem"
    error('todo')
    [x1d,custom,dxMin] = setupSEMCustom(xminmax, dx1, scheme_order, sem_iface);
    
    dt_min = dxMin/(c*cfl);
    if dt > dt_min
        error('dt > dt_min')
    end
    
    sourceFactor1 = sourceScalingSEM(scheme_order);
    
    domain1 = Domain1D(x1d, dx1, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fmax,sourceFactor1,domain1);
    s1 = Simulation1D(SolverType.SEM, domain1, src1, custom);
    
    p = runSingleDomainSolver(iter, s1);    
    p = p;
    
    if write_data
        save(sprintf('%s/%s_L%0.2f_dx%0.3f_dt%0.6f_p%i.mat', base_path, conv_test, L, dx1, dt, scheme_order), ...\
            'p', 'x1d', 'dx', 'dt', 'tmax', 'l', 'scheme_order');
    end
elseif conv_test == "fourier_sem"  
    xminmax_1 = [0,l_d1];
    xminmax_2 = [0,l_d2];
    
    source_partition = "LEFT";
    [p1,p2,s1,s2,iter] = run_FOURIER_SEM(tmax, xminmax_1, xminmax_2, dx1, dx2, ...\
        fmax, source_partition, x0_pos, c, rho, scheme_order, interface_order, iface);
    
    dt = s1.domain.dt;
    %norm_factor = 1.0/max([p1(iter,:), p2(iter,:)]);
    
    L = s1.domain.l + s2.domain.l;
    
    x1d_1 = s1.domain.x1d;
    x1d_2 = s2.domain.x1d;
    
    dx = dx2;
    
    if write_data
        save(sprintf('%s/%s_L%0.2f_dx%0.3f_dt%0.6f_p%i_ip%i.mat', ...\
            base_path, conv_test, L, dx2, dt, scheme_order, interface_order), ...\
            'p1', 'p2', 'x1d_1', 'x1d_2', 'dx', 'dt', 'tmax', 'L', 'scheme_order', 'interface_order');
    end
else
    error('test case not supported')
end

% calc RMS

xminmax = [0,L];
[p_ref, x1d_ref] = calcAnalytical(xminmax, dt, dx1, fmax, x0_pos, L, c, rho, tmax);

% indx = int32(tmax/dt);
% x1d_full1 = [x1d_1(1:end-1), l_d1+x1d_2];    
% 
% p_ref = spline(x1d_ref, p_ref, x1d_full1);
% 
% p1_ = p1(indx,:); p2_ = p2(indx,:);
% err = abs(p_ref(indx,:) - [p1_(1:end-1),p2_])./p_ref(indx,:);

%% PLOT
if ~do_plot
    return
end

if conv_test == "fourier_sem"
    h = figure(1);
    set(h,'Position',[10 10 1600 800])
    tiledlayout(2,1)
    ax1 = nexttile;
    ax2 = nexttile;
else
    h = figure(1);
    set(h,'Position',[10 10 1600 800])
end

for n=1:iter    
    t = n*dt;
    if mod(n-1,5) == 0 || n == iter
        if conv_test == "fourier_sem"
            title_str = sprintf('%s | %s \n ppw ratio %i : %i          dt ratio 1 : %i \n t=%f', .../
                s1.solver_type, s2.solver_type, ppw_fourier, ppw_sem, s1.domain.dt/s2.domain.dt, t);
            plotSemiLogCoupled(x1d_1, l_d1+x1d_2, x1d_ref, p1(n,:), p2(n,:), p_ref(n,:), title_str, ax1, ax2)
            drawnow
        else
            title_str = sprintf('t = %0.4f, n = %i', t, n);
            p_n = p(n,:);

            if mod(n-1, 5) == 0 || n == iter
                plot(x1d,p_n,'-o')        
                hold on
                hold off
                xlabel('x')
                ylabel('p')
                ylim([-1,1])
                xline(l_d1,'--')
                title(title_str)
                drawnow
            end
        end

        if write_gif
            filename = sprintf('%s_L%0.2f_dx%0.3f_dt%0.6f.gif', conv_test, L, dx1, dt);
            writeGif(h, sprintf('%s/%s',base_path, filename), n==1);
        end
    end
end


function [p, x1d] = calcAnalytical(xminmax, dt, dx1, fmax, x0_pos, L, c, rho, tmax)
    iter = tmax/dt;
    if rem(iter,1) ~=0
        error('ERROR: tmax not a multiple of dt')
    end
    
    nmodes = L/(2*dx1);
    x1d = xminmax(1):dx1:xminmax(end);
    
    sourceFactor = sourceScalingFourier();
    
    domain1 = Domain1D(x1d, dx1, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fmax,sourceFactor,domain1,x0_pos);
    s1 = Simulation1D(SolverType.FOURIER, domain1, src1, CustomFourier(nmodes));
    
    p = runSingleDomainSolver(iter, s1);
    p = p*2;
end