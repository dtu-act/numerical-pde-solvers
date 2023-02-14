clear all
close all
addpath('../shared')

base_path = '/Users/nikolasborrel/data/deeponet/input_1D_2D/matlab';
if ~exist(base_path, 'dir')
    status = mkdir(base_path);

    if status ~= 1
        error("Error writing dir. Error code: %i",status)
    end
end

write_data = true;
plot_impedance_fitting = true;

%% SIMULATION PARAMETERS
source_type = "gaussian"; % "gaussian" | "grf"
boundary_type = "freq_indep"; % "neumann", "dirichlet", "freq_indep" | "freq_dep"

tmax = 18;
xminmax = [0,5];

c_phys = 343.0; % m/s speed of sound
fmax_phys = 1000;
rho = 1.2;

src_padding = 0.6; % distance to the boundaries

% resolutions
ppw = 8; % points per wavelength for SEM simulation
ppw_u = 8; % ppw for branch net u (should be set low for GRF)
ppw_srcs = 5; % resolution for the source samples distributions

% resolution for output data
ppw_x_out = 4;
ppw_t_out = 4;

Porder = 4;

c = 1;
fmax = fmax_phys/c_phys; % Hz
sigma0 = c/(pi*fmax/2);

L = xminmax(2)-xminmax(1);
num_of_nodes_u = ceil(L/(c_phys/(fmax_phys*ppw_u)));
num_srcs = ceil(L/(c_phys/(fmax_phys*ppw_srcs)));

if boundary_type == "freq_dep"
    %sigma = 8000/c_phys; % Flow resistivity
    %d = 0.1; % Thickness of porous material inside cavity
    sigma = 8000/c_phys;
    dmat = 0.1;
    filter_order = 4; % Order of approximation on rational function for BC's
    CFL = 0.1;
elseif boundary_type == "freq_indep"
    %xi = 5.83; % specific acoustic impedance
    xi = 17.98; % specific acoustic impedance
    CFL = 1.0;
else
    CFL = 1.0;
end

%% SETUP SIMULATION
[vx,etov,conn,rs,x1d] = setup.setupMesh1D(c,fmax,Porder,ppw,xminmax);
dxMin = min(x1d(2:end) - x1d(1:end-1)); % assume sorted
dt = CFL*dxMin/c; % courant condition
timesteps = ceil(tmax/dt);

if source_type == "gaussian"
    x0_srcs = linspace(xminmax(1)+src_padding,xminmax(2)-src_padding,num_srcs);
    %x0_srcs = [-0.3,-0.15,0.0,0.15,0.3];
    x1d_u = linspace(x1d(1), x1d(end), num_of_nodes_u);
    [p0_ics, x_ic] = generateGaussianICs(x1d_u, x0_srcs, c, sigma0);
elseif source_type == "grf"
    x0_srcs = [];
    [p0_ics, x_ic] = generateGRFs(xminmax, fmax, c, [0.3], num_srcs, num_of_nodes_u);
else
    error('ic not supported')
end

num_accs = 4;
accs_srcs = zeros(size(p0_ics,1),timesteps,num_accs*2); % left and right accs

figure(1)
plotICs(x_ic, p0_ics)

p_hat_srcs = zeros(size(p0_ics,1), timesteps, length(x1d));
v_hat_srcs = p_hat_srcs;

dtAccs = @(accs,p) 0;
if boundary_type == "neumann"
    impedance_data.type = boundary_type;
    bc_update_fn = @(rhs_p,p,accs) rhs_p; % no updates
elseif boundary_type == "dirichlet"
    impedance_data.type = boundary_type;
    bc_update_fn = @(rhs_p,p,accs) 0;
elseif boundary_type == "freq_indep"
    impedance_data.type = boundary_type;
    impedance_data.Z = xi*rho*c;
    Z = impedance_data.Z;    
    bc_update_fn = @(rhs_p,p,accs) -rho*c^2*p/Z + rhs_p;    
elseif boundary_type == "freq_dep"    
    f_range = (50:1:5000)/c_phys; % Range where boundary conditions are defined
    [impedance_data, Y, fit, fspan] = impedance_bcs.fitImpedanceBoundaryData(rho, c, f_range, sigma, filter_order, dmat);
    impedance_data.type = boundary_type;

    data = impedance_data;
    if plot_impedance_fitting
        impedance_bcs.plotImpedanceFitting(fspan, c_phys, Y, fit)
    end
    
    dtAccs = dtAccsSetup(data.lambdas,data.alpha,data.beta);
        
    % velocity at the boundaries
    v_fn = @(p,accs) data.Yinf*p + data.A(1)*accs(1) + data.A(2)*accs(2) + 2*(data.B*accs(3)+data.C*accs(4));    
    bc_update_fn = @(rhs_p,p,accs) -rho*c^2*v_fn(p,accs) + rhs_p;   
else
    error('boundary type not supported')
end

[A, Sx] = setup.assembleWaveEquationCoupled1D(vx,etov,conn,rs);

%% RUN SIMULATION
tic
for i=1:size(p0_ics,1)
    p0_fn = @(x) spline(x_ic, p0_ics(i,:), x);
        
    [p_hat, v_hat, laccs, raccs] = solvers.solveCoupledWaveEquation1D(...
        tmax, c, rho, p0_fn(x1d), @(t) 0, A, Sx, timesteps, bc_update_fn, dtAccs);
    
    % normalize pressure
    pmax = max(abs(p_hat), [], 'all');
    
    p_hat_srcs(i,:,:) = p_hat/pmax;
    v_hat_srcs(i,:,:) = v_hat/pmax;
    
    accs_srcs(i,:,:,:,:) = [laccs,raccs];
end

fprintf('Simulation time total: %0.2f\n', toc)
fprintf('Dof: %d\nn: %d\ndt_sim: %1.8f\n\n',size(A,1),timesteps,dt);

tsteps = linspace(0,tmax,size(p_hat_srcs, 2));

if write_data
    if boundary_type == "freq_indep"
        path_file = sprintf('%s/%s_1D_%0.2fHz_sigma%0.1f_c%i_xi%0.2f_%s%i_T%i.h5',...
            base_path,boundary_type,fmax*c_phys,sigma0,c,xi,source_type,size(p0_ics,1), tmax);
    elseif boundary_type == "freq_dep"
        path_file = sprintf('%s/%s_1D_%0.2fHz_sigma%0.1f_c%i_d%0.2f_%s%i_T%i.h5', ...
            base_path,boundary_type,fmax*c_phys,sigma0,c,dmat,source_type,size(p0_ics,1), tmax); 
    else
        path_file = sprintf('%s/%s_1D_%0.2fHz_sigma%0.1f_c%i_%s%i_T%i.h5', ...
            base_path,boundary_type,fmax*c_phys,sigma0,c,source_type,size(p0_ics,1), tmax);
    end
    
    if ppw_x_out > 1        
        xprune = round(ppw/ppw_x_out);
        x1d = x1d(1:xprune:end);
        p_hat_srcs = p_hat_srcs(:,:,1:xprune:end);
    end

    if ppw_t_out > 1
        tprune = round((ppw/ppw_t_out)/CFL)*2;

        tsteps = tsteps(1:tprune:end);
        p_hat_srcs = p_hat_srcs(:,1:tprune:end,:);
    end

    % transpose for .h5 format to be correct
    p_hat_srcs_ = permute(p_hat_srcs,[3,2,1]);
    v_hat_srcs_ = permute(v_hat_srcs,[3,2,1]);
    accs_srcs_ = permute(accs_srcs,[4,3,2,1]);

    write.writeHDF5(x1d,tsteps,p_hat_srcs_,conn,x0_srcs,accs_srcs_,...
        c,c_phys,rho,sigma0,fmax,ppw,boundary_type,dxMin,xminmax,path_file)
    
    if boundary_type == "freq_dep"
        path_json = sprintf('%s/freq_dep_params_d%0.4f.json',base_path,dmat);
        [filepath,name,ext] = fileparts(path_file);
        write.writeImpedanceParams(impedance_data,accs_srcs_(1,:,:),c,boundary_type,name,path_json)
    end
end

srcs_indx = 1;
p_hat = squeeze(p_hat_srcs(srcs_indx,:,:));

%% PLOTTING
figure(2)
contourf(x1d, tsteps, p_hat )
colorbar
xlabel('x [m]')
ylabel('time [sec]')   

figure(3)
for i = 1:length(tsteps)
    if mod(i-1,1) == 0
        t = tsteps(i);
        plot(x1d,p_hat(i,:))
        legend('Numerical')
        xlabel('x [m]')
        ylabel('p')
        ylim([-1,1]);          
        drawnow
    end
end

function [p0_ics, x_ic] = generateGaussianICs(x1d, x0_srcs, c, sigma0)
    p0_ics = zeros(length(x0_srcs), length(x1d));
    for i = 1:length(x0_srcs)
        source = sources.gaussianSourceIC1D(c,sigma0); % precalc sigma0 to align with previous results
        p0_ics(i,:) = source(x1d,0,x0_srcs(i));
    end
    x_ic = x1d;
end

function [p0_ics, x_ic] = generateGRFs(xminmax, fmax, c, l_0, num_of_samples, num_of_nodes)
    sigma0_window = c/(pi*fmax/2);
    [p0_ics, x_ic] = ics.generateGRFs1D(xminmax, [1.0], l_0, num_of_nodes, num_of_samples, sigma0_window);
end

function plotICs(x_ic, p0_ics)   
    for i = 1:size(p0_ics,1)    
        plot(x_ic, p0_ics(i,:))
        hold on
    end
    hold off

    xlabel('x')
    ylabel('Pressure [Pa]')
    set(gca,'fontsize',16)
end

function dtAccs_ = dtAccsSetup(lambda,alpha,beta)
    function accs = dtAccs(accs, p)
            % Boundary ADE - eq. (7) from Finnur 2019
            dphi1 = -lambda(1)*accs(1) + p;
            dphi2 = -lambda(2)*accs(2) + p;
            dpsi1 = -alpha(1)*accs(3) - beta(1)*accs(4) + p;
            dpsi2 = -alpha(1)*accs(4) + beta(1)*accs(3);
            
            accs = [dphi1;dphi2;dpsi1;dpsi2];
    end

    dtAccs_ = @dtAccs;
end