clear all
close all

do_animation = true;
write_gif = do_animation && false;

write_data = false;
write_fig = false;

base_path = '~/data/fm-sem';

%% SETUP SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%tmax = 0.0135;
tmax = 0.2;

coupling_methods = "SEM"; % SEM | FOURIER
x0_pos = 2.5;

% SEM SETUP
P_order = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xminmax = [0,13]; % domain length
fmax = 1000;     % [Hz]
c = 343;         % m/sec speed
rho = 1.2;
cfl = sqrt(3);

%% SETUP and SOLVE
switch coupling_methods    
    case "SEM"
        solver_type = SolverType.SEM;
        ppw = 8;
        [dx] = calcSpatial(xminmax(2),fmax,ppw,c);        
        
        % adjust SEM length for integer elements
        xminmax(2) = ceil(xminmax(2)/(dx*P_order))*(dx*P_order);
        fprintf('Adjusted SEM domain length: %f\n', xminmax(2))
        
        [x1d,custom,dxMin] = setupSEMCustom(xminmax, dx, P_order, {});    
        dt = dxMin/(c*cfl);
        
        iter = round(tmax/dt);

        sourceFactor = sourceScalingSEM(P_order);

        bound_type = BoundaryType.Neumann;
        domain = Domain1D(x1d, dx, dt, c, rho, xminmax, bound_type);
        src = defaultSource(fmax,sourceFactor,domain,x0_pos);
        s1 = Simulation1D(SolverType.SEM, domain, src, custom);
    case "FOURIER"
        solver_type = SolverType.FOURIER;
        ppw = 4;  
        [dx,nmodes] = calcSpatial(xminmax(2),fmax,ppw,c);
        custom = CustomFourier(nmodes);
        
        x1d_1 = xminmax_1(1):dx:xminmax_1(2);
        
        dt = dx/(c*cfl);
        iter = round(tmax/dt);
        
        sourceFactor = sourceScalingFourier(P_order);
        
        bound_type = BoundaryType.Neumann;
        domain = Domain1D(x1d_1, dx, dt, c, rho, xminmax, bound_type);
        src = defaultSource(fmax,sourceFactor,domain);
        s1 = Simulation1D(SolverType.FOURIER, domain, src, custom);
    otherwise
        error('method not supported');
end

p = s1.p_history;
F = s1.src.F;
p_all = zeros(iter,length(s1.p_current));

tstart_sim_tot = tic;
for n=1:iter
    [p,~] = solverWE1D_o(n,solver_type,p,xminmax,c,dt,custom,F,bound_type);            
    p_all(n+1,:) = p(:,1);
    
%     if mod(n-1,100) == 0
%         fprintf('Calculating n=%i ...\n',n-1)
%     end
end
sim_tot = toc(tstart_sim_tot);

x1d = s1.domain.x1d;

% normalize
tnorm = 0.007;
n_norm = round(tnorm/s1.domain.dt);
norm_factor = max(p_all(n_norm,:));

p_all = p_all/norm_factor;

% ref
tstart_ref_tot = tic;
[p_ref, s_ref] = calcReference(xminmax, ppw, dt, fmax, x0_pos, c, rho, iter, tnorm);
ref_tot = toc(tstart_ref_tot);
x1d_ref = s_ref.domain.x1d;

if write_data
    L = x1d(2) - x1d(1);
    save(sprintf('%s/%s_ppw%i_%i_src_%s_fmax_%i.mat', ...\
        base_path, coupling_methods, ppw, ppw2, source_partition, fmax), ...\
        'p1', 'p2', 'p_ref', 's1', 's2', 's_ref', 'ppw1', 'ppw2', 'tmax', 'L', 'interface_order', 'fmax');
end

fprintf('SIM TOTAL: %0.3f \n', sim_tot)
fprintf('REF TOTAL: %0.3f \n', ref_tot)
fprintf('FACTOR SIM/REF: %0.3f \n', sim_tot/ref_tot)

%% PLOT
title_str = sprintf('%s \n ppw %i   dt %i', s1.solver_type, ppw, s1.domain.dt);
legend1 = sprintf('%s',s1.solver_type);

if do_animation
    h = figure(2);
    set(h,'Position',[10 10 1600 800])
    ax1 = gca;
    
    for n=1:iter
        if mod(n-1,5) == 0
            t = n*s1.domain.dt;
            fprintf('n=%i, t=%0.4f\n',n,t)            

            plot(ax1,x1d,p_all(n,:),'o','LineWidth', 2)
            hold(ax1, 'on')
            plot(ax1,x1d_ref,p_ref(n,:), '--')
            hold(ax1, 'off')
            legend(ax1, legend1, 'Reference')
    
            xlabel(ax1,'x')
            ylabel(ax1,'p')
            title(ax1, title_str)
            ylim(ax1,[-1,1])
            xlim(ax1,[x1d(1),x1d(end)])    
            set(ax1,'FontSize',20)

            drawnow

            if write_gif
                path = sprintf('%s/%s_ppw%i_%i_fmax_%i',base_path,coupling_methods,ppw,ppw2,fmax);
                writeGif(h, path, n==1);
            end
        end
    end    
end 