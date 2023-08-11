clear all, clc, 
close all
addpath('../shared')
addpath('../SEMSolvers')

grid_res = 0.1;

couple_solutions = true;

% === AcousticParameters === %
f = 300;
acoustic_params = models.AcousticParameters(f,343.0,1.225);
sim_params = wbm.SimulationParameters2D(4);

Nw_greens = 200;

bound_cond = models.BoundaryCondition.Velocity;

% === DOMAIN 1 === %
lx1 = 4;
ly1 = 2;
offset1 = [0 0];
gcoord1 = [0.0 0.0; lx1*1.0 0.0; lx1*1.0 ly1*1.0; 0.0 ly1*1.0] + offset1; % coordinates at bound
origin1 = min(gcoord1);

bound_nodes1        = [1 2; 3 4; 4 1]; % boundary definition
bound_nodes1_single = [1 2; 2 3; 3 4; 4 1]; % boundary definition
interface_nodes1 = [2 3]; % boundary definition

r1 = [lx1/3,ly1*1/4]+origin1;
% point source
point_source1 = wbm.PointSource(r1,1);

source_line1 = wbm.SourceModel(models.SourceType.Line, "N/A", [0.0; 0.0; 1.0], [0 0; 0 0; 1 0]);
source_line1_single = wbm.SourceModel(models.SourceType.Line, "N/A", [0.0; 0.0; 0.0; 1.0], [0 0; 0 0; 0 0; 1 0]);

% === DOMAIN 2 === %
lx2 = 2;
ly2 = 6;
offset2 = [lx1 0];
gcoord2 = [0.0 0.0; lx2*1.0 0.0; lx2*1.0 ly2*1.0; 0.0 ly2*1.0; 0.0 ly1*1.0] + offset2; % coordinates at bound
origin2 = min(gcoord2);

bound_nodes2 = [1 2; 2 3; 3 4; 4 5]; % boundary definition
bound_nodes2_single = [1 2; 2 3; 3 4; 4 1]; % boundary definition
interface_nodes2 = [5 1]; % boundary definition

% point source (for testing single domains only)
r2 = [lx2/3,ly2/3] + origin2;
point_source2 = wbm.PointSource(r2,1);

% === BUILD domains === %
domain_single1 = wbm.Domain2D(lx1, ly1, gcoord1, bound_nodes1_single, ...
    origin1, bound_cond, point_source1);

domain_single2 = wbm.Domain2D(lx2, ly2, gcoord2, bound_nodes2_single, ...
    origin2, bound_cond, point_source2);

domain1 = wbm.DomainCoupling2D(lx1, ly1, gcoord1, bound_nodes1, ...
    interface_nodes1, origin1, bound_cond, point_source1);

domain2 = wbm.DomainCoupling2D(lx2, ly2, gcoord2, bound_nodes2, ...
    interface_nodes2, origin2, bound_cond);

% === CALCULATE AND PLOT === %

if couple_solutions
    domains = {domain1, domain2};
    plotting.geometryPlot(domains,3)
    
    disp("solveHelmholtz2D....")
    [P, X, Y] = wbm.solveHelmholtz2D(domains, acoustic_params, sim_params, grid_res);
    
    plotting.plotHelmholtzDomains(P,X,Y);    
else
    domains = {domain_single1, domain_single2};
    plotting.GeometryPlot(domains,3)
    
    disp("solveHelmholtz2D....")
    [P, X, Y] = wbm.solveHelmholtz2D(domains, acoustic_params, sim_params, grid_res);

    Pref = greens.greensSpatialDomains2D(X, Y, domains, acoustic_params.k);
    plotting.plotHelmholtzError(P,Pref,X,Y,'Greens','WBM');    
end