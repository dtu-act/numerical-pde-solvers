clear all
close all
addpath('../shared')

hpc_env = false;
base_path = setupEnv(hpc_env);

%% USER PARAMETERS
do_plots = true;
do_write_data = true;
plot_impedance_fitting = false;

source_type = "gaussian"; % "gaussian" | "grf"
boundary_type = "neumann"; % "neumann" | "dirichlet" | "freq_indep" | "freq_dep"

xminmax = [0.0,2.0];
yminmax = [0.0,2.0];

fixed_source_pos = true;

if fixed_source_pos 
%     x0_srcs = [[0.9,0.9];
%                [1.2,1.2];
%                [1.5,1.5];
%                [1.8,1.8];
%                [2.1,2.1]];

%     x0_srcs = [[0.7, 0.7];
%                [0.85,0.85];
%                [1.0, 1.0];
%                [1.15,1.15];
%                [1.30,1.30]];

x0_srcs = [[0.7, 0.7]];

end

tmax = 17; % normalized
fmax_phys = 1000;
c_phys = 343.0; % m/s speed of sound
rho = 1.2;

src_pad = 0.6; % distance to the boundaries

% resolutions
ppw = 4; % points per wavelength for SEM simulation
ppw_u = 2; % ppw for branch net u (should be set low for GRF)
ppw_srcs = 5; % resolution for the source samples distributions

% resolution for output data
ppw_x_out = 4;
ppw_t_out = 4;

Porder = 4;

%%%%%%%%%%%%%%%
N_accs = 4; % HARDCODED - do not change

c = 1;
fmax = fmax_phys/c_phys; % Hz
sigma0 = c/(pi*fmax/2); % 0.22 corresponds to 1000 Hz;

if ~exist(base_path, 'dir')
    status = mkdir(base_path);

    if status ~= 1
        error("Error writing dir. Error code: %i",status)
    end
end

NODETOL = 1e-8;
Np = (Porder+1)*(Porder+2)/2; % number of nodes in each element needed to support order P basis functions

if boundary_type == "freq_dep"
    sigma = 10000/c_phys; % Flow resistivity
    d = 0.03; % Thickness of porous material inside cavity
    filter_order = 6; % Order of approximation on rational function for BC's
    CFL = 0.2;
elseif boundary_type == "freq_indep"
    xi = 5.83; % specific acoustic impedance (used for freq. independent boundaries only)
    CFL = 1.0;    
else
    CFL = 1.0;
end   

%% SETUP SIMULATION
[conn,X2D,Y2D,VX,VY,etov,etoe,etof,x,y,r,s,Nk,gidx,dx] = setup.setupMesh2D(c,fmax,Porder,ppw,xminmax,yminmax,do_plots);
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
    [impedance_data, Y, fit, fspan] = impedance_bcs.fitImpedanceBoundaryData(rho, c, f_range, sigma, filter_order, d);
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
timesteps = ceil(tmax/dt);
dt = tmax/timesteps;

[X2D_u,Y2D_u] = generateRectilinearGrid(xminmax,yminmax,c,fmax,ppw_u);

if source_type == "gaussian"
    fixed_source_pos = true;

    if ~fixed_source_pos
        wavelength_min = c/fmax;
        dx_srcs = wavelength_min/ppw_srcs;

        [VXsrcs, VYsrcs] = setup.MeshGenDistMesh2D(xminmax+src_pad,yminmax-src_pad,dx_srcs,do_plots);
        x0_srcs = [VXsrcs' VYsrcs'];
    end
        
    up0_ics = generateGaussianICs2D(X2D_u,Y2D_u,x0_srcs,sigma0); % samples for branch net saved to disk only
    p0_ics = generateGaussianICs2D(X2D,Y2D,x0_srcs,sigma0); % initial condition for simulation
    
    if do_plots
        figure()
        plotICs2D(X2D,Y2D,p0_ics(1,:))
        figure()
        plotICs2D(X2D_u,Y2D_u,up0_ics(1,:))
    end    
elseif source_type == "grf"
    num_grf_samples = 2;              
    x0_srcs = [];
    % calculating GRFs is expensive, use scarse number of points
    p0_ics = ics.generateGRFsUnStructured2D([1.0], [0.3], X2D_u, Y2D_u, dx, num_grf_samples, sigma0);

    if do_plots
        figure()
        plotICs2D(x2d_u,y2d_u,p0_ics(1,:))
    end
else
    error('ic not supported')
end

% TODO accumulator data save
%laccs_srcs = zeros(size(p_ics,1),timesteps,N_accs);
%raccs_srcs = zeros(size(p_ics,1),timesteps,N_accs);    

p_all = zeros(size(p0_ics,1), length(X2D), timesteps+1);
vx_all = p_all; vy_all = p_all;

%% RUN SIMULATION
tic
%parfor i=1:size(srcs,1)
for i=1:size(p0_ics,1)
    fprintf('Calculating solution %i/%i\n', i, size(p0_ics,1))
    z0 = zeros(3*N,1);
    
    if source_type == "grf"
        % beacuse GRFs are computational expensive, we interpolate from
        % sparse data
        p_ic = griddata(X2D_u,Y2D_u,p0_ics(i,:),X2D,Y2D,'cubic');
        assert(~any(isnan(p_ic)))  
    else
        p_ic = p0_ics(i,:);
    end

    z0(2*N+1:3*N,1) = p_ic;

    [p_i,vx_i,vy_i] = solvers.solveCoupledWaveEquation2D(z0,tmax,rho,N,M2D,Sx,Sy, ...
        timesteps,bc_update_fn,accs);
        
    pmax = max(abs(p_i), [], 'all'); % normalize pressure    
    assert(pmax > 0)

    p_all(i,:,:) = p_i/pmax;
    vx_all(i,:,:) = vx_i/pmax;
    vy_all(i,:,:) = vy_i/pmax;
end
fprintf('Simulation time total: %0.2f\n', toc)
fprintf('Dof: %d\nn: %d\ndt_sim: %1.8f\n\n',size(M2D,1),timesteps,dt);

tsteps = linspace(0,tmax,size(p_all, 3));
mesh = [X2D(:) Y2D(:)];
umesh = [X2D_u(:) Y2D_u(:)];
umesh_shape = size(X2D_u);

if do_write_data    
    if boundary_type == "freq_indep"
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_xi%0.2f_%s%i.h5',...
            base_path,boundary_type,fmax*c_phys,sigma0,c,xi,source_type,size(p0_ics,1));
    elseif boundary_type == "freq_dep"
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_d%0.2f_%s%i.h5', ...
            base_path,boundary_type,fmax*c_phys,sigma0,c,d,source_type,size(p0_ics,1)); 
    else
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_%s%i.h5', ...
            base_path,boundary_type,fmax*c_phys,sigma0,c,source_type,size(p0_ics,1));
    end

    [mesh,p_out,dx_out] = write.pruneSpatial(mesh,p_all,dx,ppw,ppw_x_out);
    [tsteps,p_out,dt_out] = write.pruneTemporal(dt,fmax,tsteps,p_out,ppw_t_out);
    
    % permute for .h5 format to be correct
    p_out_perm = permute(p_out,[3,2,1]);
    up_ics_perm = permute(up_ics,[3,2,1]);  

    write.writeHDF5(mesh,umesh,umesh_shape,p_out_perm,up0_ics_perm,tsteps,conn,x0_srcs,{},...
        c,c_phys,rho,sigma0,fmax,boundary_type,dx_out,[xminmax' yminmax'],path_file)
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

function plotICs2D(X2D,Y2D,p0_ics)
    N = size(p0_ics,1);    
    for i = 1:N
        subplot(1,N,i)
        tri = delaunay(X2D(:),Y2D(:));
        trisurf(tri, X2D(:), Y2D(:), p0_ics(i,:));
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
        base_path = '/users/nborrelj/data/nborrelj/deeponet/input/matlab';
        parpool('threads')
        addpath("../shared")
        addpath("../distmesh")
    else
        %delete(gcp('nocreate'))
        %parpool('threads')
        base_path = '/Users/nikolasborrel/data/deeponet/input/matlab';    
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