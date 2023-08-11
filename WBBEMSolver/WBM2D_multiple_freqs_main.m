clear all, clc, 
close all
addpath('../shared')
addpath('../SEMSolvers')

[dir_cache, ~, dir_data] = paths.setupPaths("HPC");

Nw_greens = 150;
grid_res = 0.1;

% === AcousticParameters === %
f_min_max = [20,1200];
lx = 2;
ly = 2;

ws_max = [4,20];

tsim = 2.0;
xy0_rel = [3/5,3/5];
xy0 = [lx,ly].*xy0_rel;

NeX = lx/grid_res; NeY = ly/grid_res;

boundary_type = models.BoundaryCondition.Velocity;

% ========= %
c = 343; rho = 1.225;
aparams_f = @(k) models.AcousticParameters(k,c,rho);
bbox = [0,0; lx,ly];

gcoord = [0.0 0.0; lx*1.0 0.0; lx*1.0 ly*1.0; 0.0 ly*1.0]; % coordinates at boundaries
origin = [0, 0];

bound_nodes = [1 2; 2 3; 3 4; 4 1];

xy0 = xy0+origin;

source = wbm.PointSource(xy0,1);

% Domains
domain = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, boundary_type, source);

for w_max=ws_max
    wx_max = w_max; wy_max = w_max;
    fprintf('******************\nCalc WBM w_max=%i....\n', w_max)
    solver_f = @(k) wbm.solveHelmholtzSingle2D(domain, aparams_f(k*c/(2*pi)), ...
        wx_max, wy_max, grid_res);

    tic
    [Pfreqs_wbm, XY, xfreqs] = validation.calculateFrequencies(solver_f,c,rho,f_min_max,tsim);
    t = toc;
    
    Nwx = 2*(w_max + 1); Nwy = 2*(w_max + 1); % same in x and y    
    caching.saveSimulationData(dir_data, Pfreqs_wbm, XY, xfreqs, t, ...
        f_min_max, 1, NeX, NeY, bbox, xy0_rel, boundary_type, aparams_f(0), sprintf('WBM_Nwx%i_Nwy%i', Nwx,Nwy))
                                                
    solver_greens_f = @(k) caching.cacheGreensK(aparams_f(k*c/(2*pi)),...
                                                dir_cache,...
                                                1,NeX,NeY,XY,...
                                                xy0_rel,...
                                                bbox,...
                                                boundary_type,...                                                
                                                Nw_greens);

    tic
    [Pfreqs_ref, XY_greens, xfreqs_ref] = validation.calculateFrequencies(...
        solver_greens_f,c,rho,f_min_max,tsim);
    t = toc;

    assert(all(XY == XY_greens, 'all')); assert(all(xfreqs == xfreqs_ref, 'all'))    
    caching.saveSimulationData(dir_data, Pfreqs_ref, XY, xfreqs, t, ...
        f_min_max, 1, NeX,NeY, bbox, xy0_rel, boundary_type, aparams_f(0), 'GREENS_REF')
    
    plotting.plotTF(xfreqs,Pfreqs_wbm,Pfreqs_ref,Nwx,XY,[0.2,0.2],"WBM")
end