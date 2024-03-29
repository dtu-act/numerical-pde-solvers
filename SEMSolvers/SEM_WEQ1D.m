clear all
close all
addpath('../shared')

base_path = '/Users/nikolasborrel/data/deeponet/input_1D_2D/matlab/';
if ~exist(base_path, 'dir')
    status = mkdir(base_path);
    if status ~= 1
        error("Error writing dir. Error code: %i",status)
    end
end

write_data = true;
do_plots = true;

%% SIMULATION PARAMETERS
source_type = "gaussian"; % "gaussian" | "grf"
boundary_type = "freq_indep"; % "neumann", "dirichlet", "freq_indep" | "freq_dep"
src_fixed = true;

tmax = 18;

xminmax = [0,3];    % trunk net input dimension (non-uniform grid)
xminmax_u = [0,3];  % branch net input dimension (uniform grid)- could differ from xminmax when used for transfer learning

c_phys = 343.0;     % [m/s] speed of sound
fmax_phys = 1000;   % max frequency
rho = 1.2;          % density of air

src_padding = 0.6;  % distance to the boundaries (depends on max freq.)

% resolutions
ppw = 8;        % spatial resolutions (points per wavelength)
ppw_u = 2;      % branch net u resolution (should be set low for GRF)
ppw_srcs = 5;   % resolution for the source sample density

% resolution for output data (pruning)
ppw_x_out = 4; % spatial
ppw_t_out = 4; % temporal

Porder = 4; % basis function polynomial order

c = 1;                      % normalized speed of sound [m/s]
fmax = fmax_phys/c_phys;    % max frequency normalized [Hz]
sigma0 = c/(pi*fmax/2);     % Gaussian source width

L = xminmax(2)-xminmax(1);                              % domain size
L_u = xminmax_u(2)-xminmax_u(1);                        % source position domain size
num_of_nodes_u = ceil(L_u/(c_phys/(fmax_phys*ppw_u)));  % number of sensors
num_srcs = ceil(L/(c_phys/(fmax_phys*ppw_srcs)));       % number of sources

if boundary_type == "freq_dep"
    sigma = 8000/c_phys;    % flow resistivity
    dmat = 0.1;             % thickness of porous material inside cavity
    filter_order = 4;       % order of approximation on rational function for BC's
    CFL = 0.2;              % Courant condition
elseif boundary_type == "freq_indep"
    xi = 17.98;             % specific/normalized acoustic impedance
    CFL = 1.0;              % Courant condition
else
    CFL = 1.0;              % Courant condition
end

%% SETUP SIMULATION
[vx,etov,conn,rs,x1d,dx] = setup.setupMesh1D(c,fmax,Porder,ppw,xminmax);

dxMin = min(x1d(2:end) - x1d(1:end-1)); % min dx - assume sorted
dt = CFL*dxMin/c;                       % dt calculated from courant condition
timesteps = ceil(tmax/dt);

if source_type == "gaussian"
    if src_fixed
        x0_srcs = linspace(src_padding,L-src_padding,5);
    else
        x0_srcs = linspace(xminmax(1)+src_padding,xminmax(2)-src_padding,num_srcs);
    end    
    x1d_u = linspace(xminmax_u(1), xminmax_u(2), num_of_nodes_u);
    [up_ics, ux0_ic] = generateGaussianICs(x1d_u, x0_srcs, c, sigma0);
    [p_ics, x0_ic] = generateGaussianICs(x1d, x0_srcs, c, sigma0);
elseif source_type == "grf"
    x0_srcs = [];
    [p_ics, x0_ic] = generateGRFs(xminmax, fmax, c, [0.3], num_srcs, num_of_nodes_u);
else
    error('ic not supported')
end

N_accs = 4; % HARDCODED - number of accumulators, do not change
accs_all = zeros(size(p_ics,1),timesteps,N_accs*2); % left and right accumulators

figure()
plotICs(x0_ic, p_ics)
figure()
plotICs(ux0_ic, up_ics)

p_all = zeros(size(p_ics,1), timesteps, length(x1d));
v_all = p_all;

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
    if do_plots
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

%%%%%%%%%%%%%%%%%%%%%%
%%% RUN SIMULATION %%%
%%%%%%%%%%%%%%%%%%%%%%
tic
for i=1:size(p_ics,1)
    if source_type == "grf"
        % GRFs are expensive to calculate (not so much in 1D), so we
        % interpolate
        p_i = spline(x0_ic, up_ics(i,:), x1d);
    else
        p_i = p_ics(i,:);
    end

    [p_hat, v_hat, laccs, raccs] = solvers.solveCoupledWaveEquation1D(...
        tmax, c, rho, p_i, @(t) 0, A, Sx, timesteps, bc_update_fn, dtAccs);
    
    % normalize pressure
    pmax = max(abs(p_hat), [], 'all');
    
    p_all(i,:,:) = p_hat/pmax;
    v_all(i,:,:) = v_hat/pmax;
    
    accs_all(i,:,:,:,:) = [laccs,raccs];
end

fprintf('Simulation time total: %0.2f\n', toc)
fprintf('Dof: %d\nn: %d\ndt_sim: %1.8f\n\n',size(A,1),timesteps,dt);

% to be consistant with 2D code (rk code differs in spatial/temporal order - todo change)
p_all = permute(p_all,[1,3,2]); 
v_all = permute(v_all,[1,3,2]);

tsteps = linspace(0,tmax,size(p_all, 3));
mesh = x1d;
umesh = x1d_u;
umesh_shape = size(x1d_u);

if write_data    
    if boundary_type == "freq_indep"
        path_file = sprintf('%s/1D_%s_%0.2fHz_L%0.2f_sigma%0.1f_c%i_xi%0.2f_%s%i_T%i.h5',...
            base_path,boundary_type,fmax*c_phys,L,sigma0,c,xi,source_type,size(p_ics,1), tmax);
    elseif boundary_type == "freq_dep"
        path_file = sprintf('%s/%s_1D_%0.2fHz_L%0.2f_sigma%0.1f_c%i_d%0.2f_%s%i_T%i.h5', ...
            base_path,boundary_type,fmax*c_phys,L,sigma0,c,dmat,source_type,size(p_ics,1), tmax); 
    else
        path_file = sprintf('%s/%s_1D_%0.2fHz_L%0.2f_sigma%0.1f_c%i_%s%i_T%i.h5', ...
            base_path,boundary_type,fmax*c_phys,L,sigma0,c,source_type,size(p_ics,1), tmax);
    end

    [mesh,p_out,dx_out] = write.prune.pruneSpatial(mesh,p_all,dx,ppw,ppw_x_out);
    [tsteps,p_out,dt_out] = write.prune.pruneTemporal(dt,fmax,tsteps,p_out,ppw_t_out);
    
    % permute for .h5 format to be correct
    p_out_perm = permute(p_out,[2,3,1]);
    up_ics_perm = permute(up_ics,[2,3,1]);
    accs_all_prem = permute(accs_all,[4,3,2,1]);

    write.writeAllHDF5(mesh,umesh,umesh_shape,p_out_perm,up_ics_perm,tsteps,conn,x0_srcs',accs_all_prem,...
        c,c_phys,rho,sigma0,fmax,boundary_type,dx_out,xminmax,path_file)
end

%% PLOTTING
srcs_indx = 1;
p_i = p_out;
if length(size(p_out)) == 3
    % more than 1 source position has been calculated - extract first for
    % plotting
    p_i = squeeze(p_out(srcs_indx,:,:));
end


if do_plots
    figure()
    for i = 1:size(p_i,2)
        plot(mesh, p_i(:,i));
        %shading flat
        %axis equal
        xlabel('x')
        ylabel('p')
        title(sprintf('Simulated, t = %1.5f s',dt_out*(i-1)))
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