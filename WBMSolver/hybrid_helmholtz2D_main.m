clear all, clc, 
close all

couple_solutions = true;
Nw_greens = 200;

% === ACOUSTIC PARAMTERS === %
acoustic_params = models.AcousticParameters(300,343.0,1.225);
boundary_type = models.BoundaryCondition.Velocity;
bbox = [0,0;2,2];

x0 = (bbox(1,1)+bbox(2,1))/2; y0 = (bbox(1,2)+bbox(2,2))/2;
xy0 = [x0,y0];

source = models.PointSource(xy0,-1);
% sigma = 0.2; % -> 1000Hz % source = setup.gauss2DSetup(c,sigma,xy0);

%==== SEM SIMULATION PARAMETERS =====%
P_order = 4;
NeX = 16;
NeY = 16;

solverSEM = solvers.HelmholtzSolver2D(acoustic_params);
solverSEM.setupMeshUniform(bbox, NeX, NeY);
solverSEM.setupSolver(P_order, boundary_type, source);
[P, XY] = solverSEM.solve();

Pref = caching.cacheGreens(P_order,NeX,NeY,XY,...
    '~/data/greens',solverSEM.acoustic_params,solverSEM.source.r0,solverSEM.mesh_info.bbox,solverSEM.boundary_type,Nw_greens);

[Proi, Pref_roi] = validation.postProcessSolutions(P,Pref,XY,bbox,xy0,0.3);
plotting.plotHelmholtzError(Proi, Pref_roi, XY(:,1), XY(:,2), 'Greens', 'SEM')

return
%==== WBM SIMULATION PARAMETERS =====%
grid_res = 0.1;
sim_params = wbm.SimulationParameters2D(4);

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
%source_line1_single = wbm.SourceModel(models.SourceType.Line, "N/A", [0.0; 0.0; 0.0; 1.0], [0 0; 0 0; 0 0; 1 0]);

% === BUILD domains === %
domain1 = wbm.DomainCoupling2D(lx1, ly1, gcoord1, bound_nodes1, ...
    interface_nodes1, origin1, bound_cond, point_source1);

domain2 = wbm.DomainCoupling2D(lx2, ly2, gcoord2, bound_nodes2, ...
    interface_nodes2, origin2, bound_cond);

% === CALCULATE AND PLOT === %

if couple_solutions
    domains = {domain1, domain2};
    disp("solveHelmholtz2D....")
    [P, X, Y] = wbm.solveHelmholtz2D(domains, acoustic_params, sim_params, grid_res);
    
    plotting.plotHelmholtzDomains(P,X,Y);    
else
    % TESTING ONLY
    domain_single1 = wbm.Domain2D(lx1, ly1, gcoord1, bound_nodes1_single, ...
        origin1, bound_cond, point_source1);

    domain_single2 = wbm.Domain2D(lx2, ly2, gcoord2, bound_nodes2_single, ...
        origin2, bound_cond, point_source2);

    domains = {domain_single1, domain_single2};
    plotting.GeometryPlot(domains,3)
    
    disp("solveHelmholtz2D....")
    [P, X, Y] = wbm.solveHelmholtz2D(domains, acoustic_params, sim_params, grid_res);

    Pref = greens.greensSpatialDomains2D(X, Y, domains, acoustic_params.k);
    plotting.plotHelmholtzError(P,Pref,X,Y,'Greens','WBM');
end