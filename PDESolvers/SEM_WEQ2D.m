clear all
close all

hpc_env = true;

%% USER PARAMETERS

sim_id = 'cube_freq_indep_ppw4_noprune';
do_plot = true & ~hpc_env;
do_write_data = true;
plot_impedance_fitting = false;

rand_src_pos = true;

source_type = "gaussian"; % "gaussian" | "grf"
boundary_type = "freq_indep"; % "neumann" | "dirichlet" | "freq_indep" | "freq_dep"
geometry = "CUBE"; % "LSHAPE" | "CUBE"
geometry_src = "CUBE"; % "LSHAPE" | "CUBE"
xminmax = [0.0,3.0];
yminmax = [0.0,3.0];

tmax = 17; % normalized
fmax_phys = 1000;
c_phys = 343.0; % m/s speed of sound
rho = 1.2;

% resolutions
ppw = 4; % points per wavelength for SEM simulation
ppw_u = 4; % ppw for branch net u (should be set low for GRF)
ppw_srcs = 5; % resolution for the source samples distributions

% resolution for output data
ppw_x_out = -1;
ppw_t_out = -1;

Porder = 4;

%%%%%%%%%%%%%%%
base_path = setupEnv(hpc_env);

filepath_out = sprintf('%s/%s.h5', base_path, sim_id);
xl = xminmax(1); xu = xminmax(2);
yl = yminmax(1); yu = yminmax(2);
src_pad = 0.7; % distance to the boundaries

if geometry == "LSHAPE"
    xslice = xminmax(2)/2;
    yslice = yminmax(2)/2;
    pv = [xl,yl; xl,yu; xl+xslice,yu; xl+xslice,yl+yslice; xu,yl+yslice; xu,yl; xl,yl];
elseif geometry == "CUBE"
    pv = [xl,yl;xu,yl;xl,yu;xu,yu];
else 
    error("geometry unknown")
end

if rand_src_pos    
    if geometry_src == "LSHAPE"    
        %% inner srcs
        xslice = xminmax(2)/2-src_pad;
        yslice = yminmax(2)/2-src_pad;
    
        pv_srcs = [xl+src_pad,yl+src_pad; xl+src_pad,yu-src_pad; 
                   xl+xslice,yu-src_pad; xl+xslice,yl+yslice; 
                   xu-src_pad,yl+yslice; xu-src_pad,yl+src_pad; xl+src_pad,yl+src_pad];   
    elseif geometry_src == "CUBE"
        pv_srcs = [xl+src_pad,yl+src_pad;xu-src_pad,yl+src_pad;xl+src_pad,yu-src_pad;xu-src_pad,yu-src_pad];
    else 
        error("geometry src unknown")
    end
else
    if geometry_src == "LSHAPE"    
        x0_srcs = [[0.9,0.75];[1.2,0.75];[1.5,0.75];[1.8,0.75];[2.1,0.75]]; % L-shape
    elseif geometry_src == "CUBE"
        %x0_srcs = [[0.7, 0.7];[0.85,0.85];[1.0, 1.0];[1.15,1.15];[1.30,1.30]]; % cube 2x2
        x0_srcs = [[0.9,0.9];[1.2,1.2];[1.5,1.5];[1.8,1.8];[2.1,2.1]]; % cube 3x3
        %x0_srcs = [[0.9,0.9];[1.2,1.2];]; % cube 3x3
    else 
        error("geometry unknown")
    end
end

N_accs = 4; % HARDCODED - do not change, only 4 accumulators are computed

c = 1;
fmax = fmax_phys/c_phys; % Hz
sigma0 = c/(pi*fmax/2); % 0.22 corresponds to 1000 Hz;

NODETOL = 1e-8;
Np = (Porder+1)*(Porder+2)/2; % number of nodes in each element needed to support order P basis functions

if boundary_type == "freq_dep"
    sigma = 10000/c_phys; % Flow resistivity
    dmat = 0.025; % Thickness of porous material inside cavity
    filter_order = N_accs; % Order of approximation on rational function for BC's
    CFL = 0.1;
elseif boundary_type == "freq_indep"
    xi = 5.83; % specific acoustic impedance (used for freq. independent boundaries only)
    CFL = 1.0;    
else
    CFL = 1.0;
end   

%% SETUP SIMULATION
[VX,VY,etov,dx] = setup.meshPolygon2D(c,fmax,Porder,ppw,pv,do_plot);
[conn,X2D,Y2D,VX,VY,etov,etoe,etof,x,y,r,s,Nk,gidx] = setup.setupMesh2D(VX,VY,etov,Porder);
[M2D,Sx,Sy] = setup.assembleWaveEquationCoupled2D(Nk,Np,conn,x,y,r,s,Porder);

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

dt = setup.calculateDt(r,s,x,y,CFL,c,Porder,NODETOL);
num_timesteps = ceil(tmax/dt);
dt = tmax/num_timesteps;

if source_type == "gaussian"
    if rand_src_pos
        wavelength_min = c/fmax;
        dx_srcs = wavelength_min/ppw_srcs;       

        [VXsrcs, VYsrcs] = setup.meshPolygon2D(c,fmax,1,ppw,pv_srcs,do_plot);
        x0_srcs = [VXsrcs' VYsrcs'];
        fprintf('Calculating number of src positions: %i\n', size(x0_srcs,1))
    end
    
    [X2D_u,Y2D_u] = generateRectilinearGrid(xminmax,yminmax,c,fmax,ppw_u);

    up_ics = generateGaussianICs2D(X2D_u,Y2D_u,x0_srcs,sigma0); % samples for branch net saved to disk (only)
    p_ics = generateGaussianICs2D(X2D,Y2D,x0_srcs,sigma0);      % initial condition for simulation
    
    if do_plot
        figure()
        plotICs2D(X2D_u,Y2D_u,up_ics(1,:))
        for i=1:5 % show first 5
            figure()
            plotICs2D(X2D,Y2D,p_ics(i,:),etov)
        end
    end    
elseif source_type == "grf"
    num_grf_samples = 2;              
    x0_srcs = [];
    % calculating GRFs is expensive, use scarse number of points
    % TODO: non-rectangular domain
    [X2D_u,Y2D_u] = generateRectilinearGrid(xminmax,yminmax,c,fmax,ppw_u);
    p_ics = ics.generateGRFsUnStructured2D([1.0], [0.3], X2D_u, Y2D_u, dx, num_grf_samples, sigma0);

    if do_plot
        for i=1:5 % show first 5
            figure()
            plotICs2D(x2d_u,y2d_u,p_ics(i,:))
        end
    end
else
    error('ic not supported')
end

tsteps = linspace(0,tmax,num_timesteps+1);
assert(dt == tsteps(2)-tsteps(1))

mesh = [X2D(:) Y2D(:)];
umesh = [X2D_u(:) Y2D_u(:)];
umesh_shape = size(X2D_u);

p_all = zeros(size(p_ics,1), length(X2D), num_timesteps+1);

write.initHDF5(filepath_out,mesh,umesh,umesh_shape,p_all,up_ics,tsteps,conn,x0_srcs,...
    c,c_phys,rho,sigma0,fmax,boundary_type,dx,[xminmax' yminmax'])

%% RUN SIMULATION
tic
%parfor i=1:size(p_ics,1)
for i=1:size(p_ics,1)
    fprintf('Calculating solution %i/%i\n', i, size(p_ics,1))
    z0 = zeros(3*N,1);
    
    if source_type == "grf"
        % beacuse GRFs are computational expensive, we interpolate from
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
    
    write.appendHDF5(filepath_out,p_i/pmax,i);
    
    %p_all(i,:,:) = p_i/pmax;
end
fprintf('Simulation time total: %0.2f\n', toc)
fprintf('Dof: %d\nn: %d\ndt_sim: %1.8f\n\n',size(M2D,1),num_timesteps,dt);

if do_plot
    figure()
    % use this if pruned in spatial dimension
    %tri = delaunay(X2D(:),Y2D(:)); 
    % otherwise
    tri = conn(:,1:3); % only works if not pruned in spatial dimension  
    for i = 1:size(p_i,2)
        trisurf(tri, mesh(:,1), mesh(:,2), p_i(:,i));
        %shading flat
        %axis equal
        xlabel('x')
        ylabel('y')
        zlabel('p')
        title(sprintf('Simulated, t = %1.5f s',dt_out*(i-1)))
        zlim([-1.0 1.0])
        colorbar()
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
    N = size(p0_ics,1);    
    for i = 1:N
        subplot(1,N,i)
        if nargin == 3
            conn = delaunay(X2D(:),Y2D(:));
        end
        trisurf(conn, X2D(:), Y2D(:), p0_ics(i,:));
        hold on
    end
    hold off

    title('Pressure initial condition')
    xlabel('x')
    ylabel('y')
    zlabel('Pressure [Pa]')
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

function [base_path] = setupEnv(is_HPC)    
    if is_HPC
        disp("HPC")
        base_path = '/work3/nibor/1TB/deeponet/input_1D_2D';
        addpath("/zhome/00/4/50173/matlab/distmesh-master")
        addpath("/zhome/00/4/50173/matlab/matrix_fitting")
        addpath("../shared")
        %parpool('Threads')
        parpool('Processes')
    else
        disp("LOCAL")
        %delete(gcp('nocreate'))
        %parpool('threads')
        base_path = '/Users/nikolasborrel/data/deeponet/input_1D_2D/matlab/';
        addpath('../shared')
    end

    if ~exist(base_path, 'dir')
        status = mkdir(base_path);
    
        if status ~= 1
            error("Error writing dir. Error code: %i",status)
        end
    end

end

function [XX,YY] = generateRectilinearGrid(xminmax,yminmax,c,fmax,ppw)
    dx = c/(fmax*ppw);

    num_of_nodes_x = ceil((xminmax(2)-xminmax(1))/dx);
    num_of_nodes_y = ceil((yminmax(2)-yminmax(1))/dx);
    x1d = linspace(xminmax(1), xminmax(2), num_of_nodes_x);
    y1d = linspace(yminmax(1), yminmax(2), num_of_nodes_y);
    [XX, YY] = meshgrid(x1d, y1d); 
end