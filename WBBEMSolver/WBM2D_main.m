clear all, clc, 
close all
addpath('../shared')
addpath('../SEMSolvers')

[path_ref_dir] = paths.setupPaths("LOCAL");

%%%
Nw_greens = 100;
grid_res = 0.1;

% === AcousticParameters === %
f = 300;
acoustic_params = models.AcousticParameters(f,343.0,1.225);
boundary_type = models.BoundaryCondition.Velocity;

w_max = 20;

lx = 2;
ly = 2;
bbox = [0,0; lx,ly];

gcoord = [0.0 0.0; lx*1.0 0.0; lx*1.0 ly*1.0; 0.0 ly*1.0]; % coordinates at boundaries
origin = [0, 0];

bound_nodes = [1 2; 2 3; 3 4; 4 1];

xy0_rel = [3/5,3/5];
xy0 = [lx,ly].*xy0_rel;
xy0 = xy0+origin;

r0 = 0.4;

NeX = lx/grid_res; NeY = ly/grid_res;
%% BOUNDARY DEFINITION
source = wbm.PointSource(xy0,1);

line_vel = [0 0; 0 0; 1 0; 0 0];
line_press  = [1.0; 0.0; 0; 0];
source_line_vel  = models.SourceModel(models.SourceType.Line, line_vel);
source_line_pres = models.SourceModel(models.SourceType.Line, line_press);

% Domains
domain_point_s   = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, ...
    boundary_type, source);

domain_line_vel  = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, ...
    models.BoundaryCondition.Velocity, source_line_vel);
domain_line_pres = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, ...
    models.BoundaryCondition.Pressure, source_line_pres);

domain = domain_point_s;

wx_max = w_max; wy_max = w_max;

disp("WBM solve....")
[P, XY, wc] = wbm.solveHelmholtzSingle2D(domain, acoustic_params, wx_max, wy_max, grid_res);

fprintf('Nw = %i\n', wc.nw);

Pref = caching.cacheGreens(1,NeX,NeY,XY,...
    path_ref_dir,acoustic_params,xy0_rel,bbox,boundary_type,Nw_greens);

fd = @(p) min(-drectangle(p,0,lx,0,ly), dcircle(p,xy0(1),xy0(2),r0));
Proi = validation.extractROI(P(:),XY,fd);
Proi_ref = validation.extractROI(Pref(:),XY,fd);

l2_err = norm(Proi-Proi_ref, 2);
fprintf('L2_err = %e\n', l2_err)

plotting.plotHelmholtzError(Proi, Proi_ref, XY(:,1), XY(:,2), 'Greens', 'WBM')
plotting.plotHelmholtzError(Proi, Proi_ref, XY(:,1), XY(:,2), 'Greens', 'WBM', 45)