clear all
close all
addpath('../shared')

env = "LOCAL";

if env == "LOCAL"
    path_json = jsondecode(fileread('../scripts/simulation_setups/simulation_setup_local.json'));
elseif env == "HPC"
    addpath('~/matlab/distmesh-master') % use this running on HPC
    path_json = jsondecode(fileread('../scripts/simulation_setups/simulation_setup_hpc.json'));
else
    error('Environment not supported')
end

[path_ref_dir] = paths.createDirs(path_json);

Nw_greens = 100;

%%%% USER PARAMETERS %%%
f = 300;
c = 343; rho = 1.225;
Zs = 7400;
acoustic_params = models.AcousticParameters(f,c,rho,Zs);
boundary_type = models.BoundaryCondition.Impedance;
bound_cond_f = @(x,y) Zs;

lx = 2;
ly = 2;

Porder = 6;
dx = 0.1;

%%%%%%%%%%%
bbox = [0,0;lx,ly];


%ppw = 20;
%[dx, NeX, NeY] = utils.calcGridResolution(f,acoustic_params.c,Porder,ppw,lx,ly);

NeX = lx/dx;
NeY = ly/dx;
assert(NeX == round(NeX) && NeY == round(NeY));

xy0_rel = [3/5,3/5];
xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);

source = models.PointSource(xy0,1);
% sigma = 0.2; % -> 1000Hz % source = setup.gauss2DSetup(c,sigma,xy0);

uniform_mesh = true;

r0 = 0.4;

% SOLVE %
solverSEM = solvers.HelmholtzSolver2D(acoustic_params);
if uniform_mesh
    solverSEM.setupMeshUniform(bbox, NeX, NeY);
    
    uref = @(Porder,NeX,NeY,XY) caching.cacheGreens(Porder,NeX,NeY,XY,...
        path_ref_dir,acoustic_params,...
        xy0_rel,bbox,boundary_type,Nw_greens);
else
    hmin = 0.05; hmax = 0.1250;
    solverSEM.setupMeshNonUniform(bbox,xy0,hmin,hmax);
    if boundary_type == models.BoundaryCondition.Velocity
        uref = @(XY,P,M) greens.greensNeumannUnstructuredSpatial2D(Nw_greens,Nw_greens,...
            XY,xy0,solverSEM.acoustic_params.k);
    else
        uref = @(XY,P,M) greens.greensNeumannStructuredSpatial2D(Nw_greens,Nw_greens,...
            XY,xy0,solverSEM.acoustic_params.k);
    end
end

solverSEM.setupSolver(Porder, boundary_type, source, bound_cond_f);
[P, XY] = solverSEM.solve();

Pref = uref(Porder,NeX,NeY,XY);

fd = @(p) min(-drectangle(p,0,lx,0,ly), dcircle(p,xy0(1),xy0(2),r0));
Proi = validation.extractROI(P(:),XY,fd);
Proi_ref = validation.extractROI(Pref(:),XY,fd);
%Proi = validation.calcNormConstantMax(Proi,Proi_ref)*Proi;

l2_err = norm(Proi-Proi_ref, 2);
fprintf('L2 err: %e\n', l2_err)

plotting.plotHelmholtzError(Proi, Proi_ref, XY(:,1), XY(:,2), 'Greens','SEM', 90)
plotting.plotHelmholtzError(Proi, Proi_ref, XY(:,1), XY(:,2), 'Greens','SEM', 45)