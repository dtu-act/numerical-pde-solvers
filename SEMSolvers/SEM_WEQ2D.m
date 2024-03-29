clear all
close all

hpc_env = false;        % run locally or on HPC system (DTU)
run_parallel = false;    % run in parallel

plot_impedance_fitting = false; % plot impedance fittings
do_plot = true & ~hpc_env;      % plot mesh, initial conditions and more

%%%%%%%%%%%%%%%%%%%
% USER PARAMETERS %
%%%%%%%%%%%%%%%%%%%

%sim_id = 'Lshape3x3_freq_indep_ppw265_train_orig';
%sim_id = 'Lshape3x3_freq_indep_ppw248_40srcpos_val_orig';
sim_id = 'Lshape3x3_freq_indep_ppw_2_4_2_5srcpos_val_orig';

boundary_type = "freq_indep"; % "neumann" | "dirichlet" | "freq_indep" | "freq_dep"
geometry = "LSHAPE"; % "LSHAPE" | "RECT" | FURNISHED3x3
source_type = "gaussian"; % "gaussian" | "grf"

Porder = 4;         % basis function polynomial order

tmax = 17;          % normalized as tmax * c_phys
fmax_phys = 1000;   % max frequency
c_phys = 343.0;     % m/s speed of sound
rho = 1.2;          % air density

%%%% trunk net input dimension (non-uniform grid) %%%%
xminmax = [0.0,3.0];
yminmax = [0.0,3.0];

%%%% branch net input dimension (uniform grid) %%%%
% might differ from trunk net input when data should be used
% for transfer learning having same dimensionality as the source model
xminmax_u = [0.0,3.0];
yminmax_u = [0.0,3.0];

%%%% resolutions %%%%
ppw_u = 2;  % ppw for branch net u (should be set low for GRF)
ppw = 4;    % spatial resolution

% use fixed source position (typically 5 positions). If true, ppw_srcs will
% be ignored.
fixed_src_pos = true;
if ~fixed_src_pos
    ppw_srcs = 5; % resolution for the source sample density
end
src_pad = 0.6; % padding distance to the boundaries (depends on the frequency fmax_phys)

% randomly sample max_num_srcs out of the total number of src positions.
% Can be used when the source position partition is e.g. narrow and the
% meshing places most of the points on the boundaries - hence, sample with
% a higher number of ppw and extract a subset randomly
max_num_srcs = -1; % -1 means no sampling (use all source positions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the temporal resolution can be hardcoded for e.g. training and validation
% data sets with different spatial resolutions (the user should make sure
% that the Courant stability condition is satisfied)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt_fixed = -1; % -1 -> hardcoded dt not used, calculated from Courant instead

% Parameters used for "A SENSITIVITY ANALYSIS ON THE EFFECT 
% OF HYPERPARAMETERS IN DEEP NEURAL OPERATORS APPLIED TO SOUND
% PROPAGATION", N. Borrel-Jensen et al..

% 3x3 rectangle
% dt used for generating the training dataset 
% with spatial resolution of ppw=6
% dt_fixed = 0.018338726833462715;

% 2x2 rectangle
% dt used for generating the training dataset 
% with spatial resolution of ppw=6 
% dt_fixed = 0.017763845350052;

% 3x3 furnished rectangle
% dt used for generating the training dataset 
% with spatial resolution of ppw=6 
% dt_fixed = 0.016159695817490;

% 3x3 L-shape
% dt used for generating the training dataset 
% with spatial resolution of ppw=6
% dt_fixed = 0.018867924528302;

% 2.5x2.5 L-shape 
% dt used for generating the training dataset 
% with spatial resolution of ppw=6
% dt_fixed = 0.018952062430323;

%%% resolution for output data %%%
% used to prune the simulation data afterwards.
% -1 means no pruning (post-processing can be done using a Python script)
ppw_x_out = -1;
ppw_t_out = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_path = setupEnv(hpc_env, run_parallel);

if run_parallel
    base_path = sprintf('%s/%s', base_path, sim_id);
    if isfolder(base_path)    
        delete(base_path + "/*")
    else
        [status, msg, msgID] = mkdir(base_path);
        if status ~= 1
            error(msg)
        end
    end
end

xl = xminmax(1); xu = xminmax(2);
yl = yminmax(1); yu = yminmax(2);

% take padding into account and create a smaller geometry where the sources
% can move freely
if geometry == "LSHAPE"
    xslice = xminmax(2)/2;
    yslice = yminmax(2)/2;
    pv = [xl,yl; xl,yu; xl+xslice,yu; xl+xslice,yl+yslice; xu,yl+yslice; xu,yl; xl,yl];
    
    if fixed_src_pos
        x0_srcs = [[0.9,0.75];[1.2,0.75];[1.5,0.75];[1.8,0.75];[2.1,0.75]]; % 3x3
    else
        % add padding corresponding to source width
        xslice = xminmax(2)/2-src_pad;
        yslice = yminmax(2)/2-src_pad;    
        pv_srcs = [xl+src_pad,yl+src_pad; xl+src_pad,yu-src_pad; 
                   xl+xslice,yu-src_pad; xl+xslice,yl+yslice; 
                   xu-src_pad,yl+yslice; xu-src_pad,yl+src_pad; xl+src_pad,yl+src_pad];
    end
elseif geometry == "RECT"
    pv = [xl,yl;xu,yl;xl,yu;xu,yu];
    if fixed_src_pos
        %x0_srcs = [[0.7, 0.7];[0.85,0.85];[1.0, 1.0];[1.15,1.15];[1.30,1.30]]; % 2x2
        x0_srcs = [[0.9,0.9];[1.2,1.2];[1.5,1.5];[1.8,1.8];[2.1,2.1]]; % 3x3
    else
        % add padding corresponding to source width
        pv_srcs = [xl+src_pad,yl+src_pad;xu-src_pad,yl+src_pad;xl+src_pad,yu-src_pad;xu-src_pad,yu-src_pad];
    end
elseif geometry == "FURNISHED3x3"
    pv = [xl,yl;
          xl+1,yl;
          xl+1,yl+0.6;
          xl+2,yl+0.6;
          xl+2,yl;
          xu,yl;
          xu,yu;
          xu-0.5,yu;
          xu-0.5,yu-0.3;
          xu-0.8,yu-0.3;
          xu-0.8,yu;
          xu-1.25,yu;
          xu-1.25,yu-0.3;
          xu-1.75,yu-0.3;
          xu-1.75,yu-0.001; % adjustment: a small values is subtracted to ensure distmesh to create a straight edge
          xu-2.2,yu;
          xu-2.2,yu-0.3;
          xu-2.5,yu-0.3;
          xu-2.5,yu;
          xl,yu;
          xl,yl];
    
    if fixed_src_pos
        x0_srcs = [[1.5, 0.6];[1.5, 0.6];[1.5, 2.4];[1.5, 0.6];[1.5, 0.6]]; % 3x3 (same as 3D)
    else
        % create rectangle inside for src positions
        src_pad_l = src_pad + 0.6;
        src_pad_u = src_pad + 0.3;
        pv_srcs = [xl+src_pad,yl+src_pad_l;xu-src_pad_u,yl+src_pad_l;xl+src_pad,yu-src_pad_u;xu-src_pad,yu-src_pad_u];
    end
else
    error("geometry unknown")
end

c = 1;                      % normalized speed of sound [m/s]
fmax = fmax_phys/c_phys;    % normalized max frequency [Hz]
sigma0 = c/(pi*fmax/2);     % E.g. 0.22 m corresponds to 1000 Hz;

N_accs = 4;                 % HARDCODED - do not change, only 4 accumulators are computed
Np = (Porder+1)*(Porder+2)/2; % number of nodes in each element needed to support order P basis functions

if boundary_type == "freq_dep"
    sigma = 10000/c_phys;   % Flow resistivity
    dmat = 0.025;           % Thickness of porous material inside cavity
    filter_order = N_accs;  % Order of approximation on rational function for BC's
    CFL = 0.1;              % Courant condition
elseif boundary_type == "freq_indep"
    xi = 17.98;             % specific acoustic impedance
    CFL = 1.0;              % Courant condition
else
    CFL = 1.0;              % Courant condition
end   

%%% SETUP SIMULATION %%%
[VX,VY,etov,dx] = setup.meshPolygon2D(c,fmax,Porder,ppw,pv,do_plot);

[conn,X2D,Y2D,VX,VY,etov,etoe,etof,x,y,r,s,Nk,gidx] = setup.setupMesh2D(VX,VY,etov,Porder);

[M2D,Sx,Sy] = setup.assembleWaveEquationCoupled2D(Nk,Np,conn,x,y,r,s,Porder);

NODETOL = 1e-8;
N = gidx;
accs = [];
if boundary_type == "neumann"    
    impedance_data.type = boundary_type;
    bc_update_fn = @(p,vx,vy,accs) rho*c^2*(Sx'*vx + Sy'*vy);
elseif boundary_type == "freq_indep"
    [hsurf, bcmap, M1D] = setup.calculateBoundaries2D(conn, X2D, Y2D, xminmax, yminmax, Porder, NODETOL);    
    impedance_data.type = boundary_type;
    impedance_data.Z = xi*rho*c;
    Z = impedance_data.Z;
    bc_update_fn = setupFrequencyIndependentBCs(rho,c,bcmap,hsurf,N,Porder,Sx,Sy,M1D,Z);
elseif boundary_type == "freq_dep"
    [hsurf, bcmap, M1D] = setup.calculateBoundaries2D(conn, X2D, Y2D, xminmax, yminmax, Porder, NODETOL);
    f_range = (50:1:5000)/c_phys; % Range where boundary conditions are defined
    [impedance_data, Y, fit, fspan] = impedance_bcs.fitImpedanceBoundaryData(rho, c, f_range, sigma, filter_order, dmat);
    impedance_data.type = boundary_type;

    data = impedance_data;
    if plot_impedance_fitting
        impedance_bcs.plotImpedanceFitting(fspan, c_phys, Y, fit)
    end
    
    DOF = size(M2D,1);
    accs_dim = N_accs*DOF;
    accs = zeros(accs_dim,1);    
    dtAccs = dtAccsSetup(data.lambdas,data.alpha,data.beta,bcmap,accs_dim);
    bc_update_fn = setupFrequencyDependentBCs(rho,c,bcmap,hsurf,N,Porder,Sx,Sy,M1D,data,dtAccs);
else
    error('boundary type not supported')
end

if dt_fixed > 0
    dt = setup.calculateDt(r,s,x,y,CFL,c,Porder,NODETOL);
    if dt_fixed > dt
        error('dt_fixed is coarser than minimum resolution')
    end
    dt = dt_fixed;
    
    num_timesteps = ceil(tmax/dt_fixed);
    tmax = num_timesteps*dt;
else
    dt = setup.calculateDt(r,s,x,y,CFL,c,Porder,NODETOL);
    num_timesteps = ceil(tmax/dt);
    dt = tmax/num_timesteps;
end

if source_type == "gaussian"
    dx_src = -1;
    if ~fixed_src_pos
        [VXsrcs, VYsrcs, ~, dx_src] = setup.meshPolygon2D(c,fmax,1,ppw_srcs,pv_srcs,do_plot);
        x0_srcs = [VXsrcs' VYsrcs'];
        if max_num_srcs > 0
            % sample a subset of the source positions 
            rng(0,'twister');
            r = randi([1 size(x0_srcs,1)],1,max_num_srcs);
            x0_srcs = x0_srcs(r,:);
        end
    end
    fprintf('Calculating number of src positions: %i\n', size(x0_srcs,1))
    
    [X2D_u,Y2D_u,dx_u] = generateRectilinearGrid(xminmax_u,yminmax_u,c,fmax,ppw_u);

    up_ics = generateGaussianICs2D(X2D_u,Y2D_u,x0_srcs,sigma0); % samples for branch net saved to disk (only)
    p_ics = generateGaussianICs2D(X2D,Y2D,x0_srcs,sigma0);      % initial condition for simulation
    
    if do_plot
        plotICs2D(X2D_u,Y2D_u,up_ics(1:5,:))
        plotICs2D(X2D,Y2D,p_ics(1:5,:),etov)        
    end    
elseif source_type == "grf"
    % Calculating Gaussian Random Fields (GRFs) is expensive, use scarse number of points
    % NOTE: only rectangular domains are currently supported
    num_grf_samples = 2;
    dx_src = -1; % TODO
    x0_srcs = [];    
    [X2D_u,Y2D_u,dx_u] = generateRectilinearGrid(xminmax,yminmax,c,fmax,ppw_u);
    p_ics = ics.generateGRFsUnStructured2D([1.0], [0.3], X2D_u, Y2D_u, dx, num_grf_samples, sigma0);

    if do_plot
        plotICs2D(x2d_u,y2d_u,p_ics(1:5,:))
    end
else
    error('source_type not supported')
end

tsteps = linspace(0,tmax,num_timesteps+1);
assert(dt == tsteps(2)-tsteps(1))

mesh = [X2D(:) Y2D(:)];
umesh = [X2D_u(:) Y2D_u(:)];
umesh_shape = size(X2D_u);

if run_parallel
    filepath_out = sprintf('%s/%s_header.h5', base_path, sim_id);
    write.parallel.writeHeaderHDF5(filepath_out,mesh,umesh,umesh_shape,tsteps,conn,...
        c,c_phys,rho,sigma0,fmax,boundary_type,dx,dx_u,dx_src,[xminmax' yminmax'])
else
    filepath_out = sprintf('%s/%s.h5', base_path, sim_id);
    p_all = zeros(size(p_ics,1), length(X2D), num_timesteps+1);
    write.single.initHDF5(filepath_out,mesh,umesh,umesh_shape,p_all,up_ics,tsteps,conn,x0_srcs,...
        c,c_phys,rho,sigma0,fmax,boundary_type,dx,dx_u,dx_src,[xminmax' yminmax'])
end

%%%%%%%%%%%%%%%%%%%%%%
%%% RUN SIMULATION %%%
%%%%%%%%%%%%%%%%%%%%%%
tic
%parfor i=1:size(p_ics,1)
for i=1:size(p_ics,1) % use this for non-parallel execution
    fprintf('Calculating solution %i/%i\n', i, size(p_ics,1))    
    z0 = zeros(3*N,1);
    
    if source_type == "grf"
        % because GRFs are computational expensive, we interpolate from
        % sparse data
        p_ic = griddata(X2D_u,Y2D_u,p_ics(i,:),X2D,Y2D,'cubic');
        assert(~any(isnan(p_ic)))  
    else
        p_ic = p_ics(i,:);
    end

    z0(2*N+1:3*N,1) = p_ic;

    [p_i,~,~] = solvers.solveCoupledWaveEquation2D(z0,tmax,rho,N,M2D,Sx,Sy, ...
        num_timesteps,bc_update_fn,accs);
        
    pmax = max(abs(p_i), [], 'all'); % normalize pressure    
    assert(pmax > 0)
    
    if run_parallel
        x0 = x0_srcs(i,:);
        up_ic = up_ics(i,:)
        filepath_out = sprintf('%s/%s_%i.h5', base_path, sim_id, i);
        fprintf('filepath_out : %s\n', filepath_out )
        write.parallel.writeHDF5(filepath_out,p_i/pmax,up_ic,x0);
    else
        write.single.appendHDF5(filepath_out,p_i/pmax,i);        
    end
end
fprintf('Simulation time total: %0.2f\n', toc)
fprintf('Dof: %d\nn: %d\ndt_sim: %1.8f\n\n',size(M2D,1),num_timesteps,dt);

if do_plot & ~run_parallel
    dt = tsteps(2)-tsteps(1);
    figure()
    if ppw_x_out > -1
        % use this if pruned in spatial dimension
        tri = delaunay(X2D(:),Y2D(:)); 
    else
        % only works if not pruned in spatial dimension
        tri = conn(:,1:3);
    end
    for i = 1:80 %1:size(p_i,2)
        trisurf(tri, mesh(:,1), mesh(:,2), p_i(:,i));
        fontsize(16,"points")
        %shading flat
        %axis equal
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('Pressure [Pa]')
        %title(sprintf('Simulated, t = %1.5f s',dt*(i-1)))
        zlim([-0.4 0.4])
        %colorbar()
        caxis([-0.3 0.3])    
        %view([180 0])
        drawnow
        %pause(0.5)
    end
end

function srcs = generateGaussianICs2D(x2d,y2d,xy0_srcs,sigma0)
    srcs = zeros(size(xy0_srcs,1), length(x2d(:)));
    for i = 1:size(xy0_srcs,1)
        source = sources.gaussianSourceIC2D(sigma0);
        srcs(i,:) = source(x2d(:),y2d(:),xy0_srcs(i,:));
    end
end

function plotICs2D(X2D,Y2D,p0_ics,conn)
    figure('Position', [0 500 1500 300])
    %title('Pressure initial condition')
    N = size(p0_ics,1);
    for i = 1:N
        subplot(1,N,i)
        if nargin == 3
            conn = delaunay(X2D(:),Y2D(:));
        end
        trisurf(conn, X2D(:), Y2D(:), p0_ics(i,:));
        xlabel('x')
        ylabel('y')
        zlabel('Pressure [Pa]')
    end
    fontsize(16,"points")
    
end

function [rhs_fn] = setupFrequencyIndependentBCs(rho,c,bcmap,hsurf,N,P,Sx,Sy,M1D,Z)
    SxT = Sx';
    SyT = Sy';

    function [rhs_updated] = frequencyIndependentBCs(p,vx,vy)
        bdyvec = zeros(N,1);
        ptmp = zeros(1,P+1);
        
        % Loop through all boundary faces
        for kk = 1:size(bcmap,1)
            Mn = (hsurf(kk)/2)*M1D;
            ptmp(:) = p(bcmap(kk,:))';
                
            % Calculate boundary contribution
            res = (ptmp./Z)*Mn;
            
            % Insert boundary contribution
            for jj = 1:size(bcmap,2)
                nidx = bcmap(kk,jj);
                bdyvec(nidx) = bdyvec(nidx) + res(jj);
            end        
        end
    
        rhs_updated = -rho*c^2*bdyvec + rho*c^2*(SxT*vx + SyT*vy);
    end
    
    rhs_fn = @frequencyIndependentBCs;
end

function [dtAccs_] = dtAccsSetup(lambda,alpha,beta,bcmap,accs_dim)
    bdynodes = unique(bcmap);
    function rhs_accs = dtAccs(accs, p)
        rhs_accs = zeros(accs_dim,1);
        for jj = 1:length(bdynodes)
            nidx = bdynodes(jj); %Global node index
            accidx = (nidx-1)*4+1; % Accumulator index
            rhs_accs(accidx) = -lambda(1)*accs(accidx) + p(nidx);
            rhs_accs(accidx+1) = -lambda(2)*accs(accidx+1) + p(nidx);
            rhs_accs(accidx+2) = -alpha*accs(accidx+2) - beta*accs(accidx+3) + p(nidx);
            rhs_accs(accidx+3) = -alpha*accs(accidx+3) + beta*accs(accidx+2);
        end
    end

    dtAccs_ = @dtAccs;
end

function [rhs_fn] = setupFrequencyDependentBCs(rho,c,bcmap,hsurf,N,P,Sx,Sy,M1D,data,dtAccs)
    SxT = Sx';
    SyT = Sy';

    function [rhs_p, rhs_accs] = frequencyDependentBCs(p,vx,vy,accs)
        rhs_accs = dtAccs(accs, p);

        bdyvec = zeros(N,1);
        
        % Loop through boundary faces, calculate the velocity at each of the nodes
        % and insert into the bdyvec contribution vector
        for kk = 1:size(bcmap,1) % loop through all boundary faces
            Mn = (hsurf(kk)/2)*M1D; % Mass matrix of current boundary face

            % Calculate normal velocity at this face
            vtmp = zeros(1,P+1);
            for jj = 1:P+1
                nidx = bcmap(kk,jj); %Global node index
                accidx = (nidx-1)*4+1; % Accumulator index
                vtmp(jj) = data.Yinf*p(nidx) + data.A(1)*accs(accidx) + data.A(2)*accs(accidx+1) ...
                    + 2*(data.B*accs(accidx+2)+data.C*accs(accidx+3));
            end
                
            % Calculate boundary contribution
            res = vtmp*Mn;
            
            %Insert boundary contribution
            for jj = 1:P+1
                nidx = bcmap(kk,jj);
                bdyvec(nidx) = bdyvec(nidx) + res(jj);
            end            
        end
    
        rhs_p = -rho*c^2*bdyvec + rho*c^2*(SxT*vx + SyT*vy);
    end
    
    rhs_fn = @frequencyDependentBCs;
end

function [base_path] = setupEnv(is_HPC,run_parallel)    
    if is_HPC
        disp("HPC")
        base_path = '/work3/nibor/1TB/deeponet/input_1D_2D';
        addpath("/zhome/00/4/50173/matlab/distmesh-master")
        addpath("/zhome/00/4/50173/matlab/matrix_fitting")
        addpath("../shared")
        if run_parallel
            %parpool('Threads')
            parpool('Processes');
        end
    else
        disp("LOCAL")
        %delete(gcp('nocreate'))
        %parpool('Threads')
        %parpool('Processes');
        base_path = '/Users/nikolasborrel/data/deeponet/input_1D_2D/matlab';
        addpath('../shared')
    end

    if ~exist(base_path, 'dir')
        status = mkdir(base_path);
    
        if status ~= 1
            error("Error writing dir. Error code: %i",status)
        end
    end

end

function [XX,YY,dx] = generateRectilinearGrid(xminmax,yminmax,c,fmax,ppw)
    dx = c/(fmax*ppw);

    num_of_nodes_x = ceil((xminmax(2)-xminmax(1))/dx);
    num_of_nodes_y = ceil((yminmax(2)-yminmax(1))/dx);
    x1d = linspace(xminmax(1), xminmax(2), num_of_nodes_x);
    y1d = linspace(yminmax(1), yminmax(2), num_of_nodes_y);
    [XX, YY] = meshgrid(x1d, y1d); 
end