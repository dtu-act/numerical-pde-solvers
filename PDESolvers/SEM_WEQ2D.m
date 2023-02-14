clear all
close all
addpath('../shared')

%% USER PARAMETERS

hpc_env = false;

do_plots = true;
do_write_data = false;
plot_impedance_fitting = false;

source_type = "gaussian"; % "gaussian" | "grf"
boundary_type = "neumann"; % "neumann" | "dirichlet" | "freq_indep" | "freq_dep"

num_grf_samples = 2;

tmax = 6; % normalized
f_phys = 1000;
c_phys = 343.0; % m/s speed of sound
rho = 1.2;

Porder = 4;

xminmax = [-1.0,1.0];
yminmax = [-1.0,1.0];

ppw = 5; % points per wavelength

%%%%%%%%%%%%%%%

c = 1;
fmax = f_phys/c_phys; % Hz
sigma0 = c/(pi*fmax/2); % 0.2;

if hpc_env
    base_path = '/users/nborrelj/data/nborrelj/deeponet/input/matlab';
    parpool('threads')
    addpath("../shared")
    addpath("../distmesh")
else
    %delete(gcp('nocreate'))
    %parpool('threads')
    base_path = '/Users/nikolasborrel/data/deeponet/input/matlab';    
end

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
    d = 0.02; % Thickness of porous material inside cavity
    filter_order = 4; % Order of approximation on rational function for BC's
    CFL = 0.1;
elseif boundary_type == "freq_indep"
    xi = 5.83; % specific acoustic impedance (used for freq. independent boundaries only)
    CFL = 1.0;    
else
    CFL = 1.0;
end   

%% SETUP SIMULATION
[conn,X2D,Y2D,VX,VY,etov,etoe,etof,x,y,r,s,Nk,gidx,dx] = setup.setupMesh2D(c,fmax,Porder,ppw,xminmax,yminmax,do_plots);
[M2D,Sx,Sy] = setup.assembleWaveEquationCoupled2D(Nk,Np,conn,x,y,r,s,Porder);

%dtAccs = @(accs,p) 0;
if boundary_type == "neumann"    
    impedance_data.type = boundary_type;
    bc_update_fn = @(p,vx,vy) rho*c^2*(Sx'*vx + Sy'*vy);
elseif boundary_type == "freq_indep"
    [hsurf, bcmap, M1D] = setup.calculateBoundaries2D(conn, X2D, Y2D, xminmax, yminmax, Porder, NODETOL);
    impedance_data.type = boundary_type;
    impedance_data.Z = xi*rho*c;
    Z = impedance_data.Z;
    bc_update_fn = setupFrequencyIndependentBCs(rho,c,bcmap,hsurf,N,Porder,Sx,Sy,M1D,Z);
elseif boundary_type == "freq_dep"
    error('not implemented yet')
    f_range = (50:1:5000)/c_phys; % Range where boundary conditions are defined
    [impedance_data, Y, fit, fspan] = sem.fitImpedanceBoundaryData(rho, c, f_range, sigma, filter_order, d);
    impedance_data.type = boundary_type;

    data = impedance_data;
    if plot_impedance_fitting
        plotImpedanceFitting(fspan, c_phys, Y, fit)
    end
    
    dtAccs = dtAccsSetup(data.lambdas,data.alpha,data.beta);
        
    % velocity at the left boundary
    v_fn = @(p,accs) data.Yinf*p + data.A(1)*accs(1) + data.A(2)*accs(2) + 2*(data.B*accs(3)+data.C*accs(4));    
    bc_update_fn = @(rhs_p,p,accs) -rho*c^2*v_fn(p,accs) + rhs_p;   
else
    error('boundary type not supported')
end


N = gidx;
dt = setup.calculateDt(r,s,x,y,CFL,c,Porder,NODETOL);
timesteps = ceil(tmax/dt);
dt = tmax/timesteps;

num_of_nodes_x = ceil((xminmax(2)-xminmax(1))/dx);
num_of_nodes_y = ceil((yminmax(2)-yminmax(1))/dx);
d_eps = 1e-8;
x1d = linspace(xminmax(1)-d_eps, xminmax(2)+d_eps, num_of_nodes_x);
y1d = linspace(yminmax(1)-d_eps, yminmax(2)+d_eps, num_of_nodes_y);

if source_type == "gaussian"
    x0_srcs = [[-0.3,-0.3];
               [-0.15,-0.15];
               [0.0,0.0];
               [0.15,0.15];
               [0.3,0.3]];

    srcs = generateGaussianICs2D(X2D,Y2D,x0_srcs,sigma0);
    
    if do_plots
        figure()
        plotICs2D(X2D,Y2D,srcs(1,:))
    end
    
elseif source_type == "grf"
    [XX_u, YY_u] = meshgrid(x1d, y1d); % calculating GRFs is expensive, use scarse number of points
    x0_srcs = [];
    srcs = ics.generateGRFsUnStructured2D([1.0], [0.3], XX_u, YY_u, dx, num_grf_samples, sigma0);

    if do_plots
        figure()
        plotICs2D(x2d_u,y2d_u,srcs(1,:))
    end
else
    error('ic not supported')
end

num_accs = 4;
laccs_srcs = zeros(size(srcs,1),timesteps,num_accs);
raccs_srcs = zeros(size(srcs,1),timesteps,num_accs);    

p_srcs = zeros(size(srcs,1), length(X2D), timesteps+1);
vx_srcs = p_srcs; vy_srcs = p_srcs;

%% RUN SIMULATION
tic
%parfor i=1:size(srcs,1)
for i=1:size(srcs,1)
    fprintf('Calculating solution %i/%i\n', i, size(srcs,1))
    z0 = zeros(3*N,1);
    
    if source_type == "grf"
        ic = griddata(XX_u,YY_u,srcs(i,:),X2D,Y2D,'cubic');
    else
        ic = srcs(i,:);
    end

    assert(~any(isnan(ic)))    
    z0(2*N+1:3*N,1) = ic;

    [p_hat,vx_hat,vy_hat] = solvers.solveCoupledWaveEquation2D(z0,tmax,rho,N,M2D,Sx,Sy, ...
        timesteps,bc_update_fn);
    
    % normalize pressure
    pmax = max(abs(p_hat), [], 'all');
    
    assert(pmax > 0)

    p_srcs(i,:,:) = p_hat/pmax;
    vx_srcs(i,:,:) = vx_hat/pmax;
    vy_srcs(i,:,:) = vy_hat/pmax;
end
fprintf('Simulation time total: %0.2f\n', toc)
fprintf('Dof: %d\nn: %d\ndt_sim: %1.8f\n\n',size(M2D,1),timesteps,dt);

tsteps = linspace(0,tmax,size(p_srcs, 2));

if do_write_data    
    if boundary_type == "freq_indep"
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_xi%0.2f_%s%i.h5',...
            base_path,boundary_type,fmax*c_phys,sigma0,c,xi,source_type,size(srcs,1));
    elseif boundary_type == "freq_dep"
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_d%0.2f_%s%i.h5', ...
            base_path,boundary_type,fmax*c_phys,sigma0,c,d,source_type,size(srcs,1)); 
    else
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_%s%i.h5', ...
            base_path,boundary_type,fmax*c_phys,sigma0,c,source_type,size(srcs,1));
    end   
    
    grid = utils.zip2(X2D,Y2D);
    t = linspace(0,tmax,timesteps+1); % TODO: use same as simulation for no mistakes
    
    domain_minmax = zeros(2,2);
    domain_minmax(:,1) = xminmax;
    domain_minmax(:,2) = yminmax;
    
    p_srcs = permute(p_srcs,[2,3,1]); % permute for .h5 format to be correct
%     v_srcs = zeros(size(p_srcs,1), size(vx_srcs,2), 2, size(vx_srcs,3));
%     for i=1:size(p_srcs,1)
%         for k=1:size(vx_srcs,3)
%             v_srcs(i,:,:,k) = utils.zip2(vx_srcs(i,:,k),vy_srcs(i,:,k));
%         end
%     end
%     v_srcs = permute(v_srcs,[3,2,4,1]);

    write.writeHDF5(grid,t,p_srcs,conn,x0_srcs,{},...
        c,c_phys,rho,sigma0,fmax,ppw,boundary_type,dx,domain_minmax,path_file)
end

%% PLOTTING
srcs_indx = 1;
p_hat = squeeze(p_srcs(srcs_indx,:,:));

if do_plots
    figure()
    %tri = delaunay(X2D(:),Y2D(:));
    tri = conn(:,1:3);
    for i = 1:timesteps+1
        trisurf(tri, X2D(:),Y2D(:), p_hat(:,i));
        %shading flat
        %axis equal
        xlabel('x')
        ylabel('y')
        zlabel('p')
        title(sprintf('Simulated, t = %1.5f s',dt*i))
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