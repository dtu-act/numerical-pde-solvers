clear all
close all

[path_ref_dir, plot_path] = paths.setupPaths("HPC");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_avg_computations = 20;
sim_methods = {{"WBM"},{"SEM",1},{"SEM",2},{"SEM",4},{"SEM",6},{"SEM",8}}; % "GREENS", "WBM", "SEM"
Porders_plot = [2,3];

ref_method = "SEM"; % "GREENS" | "SEM"
Nw_greens = 4000; % used only when comparing with Greens

roi_type = "fd_fullnobounds";
bthickness = 0.01;
slice_source = 0.2;

%% === SETUP ACOUSTIC PARAMETERS === %

f = 300;
lx = 2;
ly = 2;
bbox = [0,0;lx,ly];

grid_res = 0.1;
r0 = 0.4;

xy0_rel = [3/5, 3/5];

c = 343.0;
rho = 1.225;
acoustic_params = models.AcousticParameters(f,c,rho);
boundary_type = models.BoundaryCondition.Velocity;

errs_all = cell(length(sim_methods),1);
cpu_timings_all = cell(length(sim_methods),1);
legend_info = cell(1,length(sim_methods));

for k=1:length(sim_methods)    
    sim_method_tuple = sim_methods{k};
    [xinputs,xlabelstr,plot_title,Porder] = setupSimulationParameters(sim_method_tuple);
    
    errs = zeros(1,length(xinputs));
    cpu_timings = zeros(1,length(xinputs));
    
    disp('*****Simulation info*****')
    fprintf('Number of avg. computations per method: %i\n', num_avg_computations)
    fprintf('Domain size: %0.1fx%0.1f\n', lx,ly)
    fprintf('Nw Greens: %i\n', Nw_greens)
    fprintf('s_{xy} = (%0.2f,%0.2f)\n\n', lx*xy0_rel(1), ly*xy0_rel(2))

    if sim_method_tuple{1} == "WBM" || sim_method_tuple{1} == "GREENS"
        xy0 = [lx,ly].*xy0_rel;

        gcoord = [0.0 0.0; lx*1.0 0.0; lx*1.0 ly*1.0; 0.0 ly*1.0]; % coordinates at boundaries
        origin = [0,0];
        bound_nodes = [1 2; 2 3; 3 4; 4 1];

        % point source
        source = wbm.PointSource(xy0,1);
        domain = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, boundary_type, source);

        %% validation plot before starting simulation %%
        fd = validation.roiFunctions(roi_type,bthickness,slice_source,bbox,xy0,r0);
        [Pval, XY] = wbm.solveHelmholtzSingle2D(domain, acoustic_params, 20, 20, grid_res);
        Proi = validation.extractROI(Pval,XY,fd);

        figure(1)
        plotting.plotHelmholtz(abs(Pval(:)), XY(:,1), XY(:,2), sim_method_tuple{1}, [1,2,1], 45)
        plotting.plotHelmholtz(abs(Proi(:)), XY(:,1), XY(:,2), sprintf('%s ROI',sim_method_tuple{1}), [1,2,2], 45)
    else            
        solverSEM = solvers.HelmholtzSolver2D(acoustic_params);
    end
        
    for j=1:num_avg_computations
        for i=1:length(xinputs)
            xinput = xinputs(i);

            tic
            switch sim_method_tuple{1}
                case "GREENS"
                    NeX = lx/grid_res; NeY = ly/grid_res;
                    Porder_ref = 6;
                    legend_info{k} = 'GREENS';                    

                    fprintf('*** Compute Greens for Nw=%i...\n', xinput)
                    assert(lx == ly)
                    P = greens.greensNeumannUnstructuredSpatial2D(xinput,xinput,XY,xy0,acoustic_params.k);
                case "WBM"
                    NeX = lx/grid_res; NeY = ly/grid_res;
                    Porder_ref = 6;
                    legend_info{k} = 'WBM';

                    fprintf('*** Compute WBM for w=%i...\n', xinput)
                    [P, XY, wavecount] = wbm.solveHelmholtzSingle2D(domain, acoustic_params, xinput, xinput, grid_res);
                case "SEM"
                    NeX = xinput; NeY = xinput;
                    max_num_nodes = (lx/grid_res*6+1)*(lx/grid_res*6+1);
                    Porder_ref = findMaxPolynomialOrder(NeX,NeY,max_num_nodes);
                    legend_info{k} = sprintf('SEM P=%i', Porder);                    

                    fprintf('*** Compute SEM for P=%i, Ne=%i...\n', Porder, xinput)

                    % point source
                    xy0 = meshing.calculateSourcePosition(xinput,xinput,bbox,xy0_rel);
                    source = models.PointSource(xy0,1);

                    solverSEM.setupMeshUniform(bbox, xinput, xinput);
                    solverSEM.setupSolver(Porder, boundary_type, source);                
                    [P, XY] = solverSEM.solve();

                    if ref_method == "SEM"  
                        nodeindxs = unique(solverSEM.conn(:,1:3));
                        P = P(nodeindxs,:);
                        XY = XY(nodeindxs,:);
                    end
                otherwise
                    error('Simulation method not supported')
            end

            cpu_timings(i) = cpu_timings(i) + toc;

            switch ref_method
                case "SEM"            
                    [Pref, XYref] = caching.cacheSEM(Porder_ref,NeX,NeY,...
                        path_ref_dir,...
                        acoustic_params,...
                        xy0_rel,...
                        bbox,...
                        boundary_type,...
                        true);
                case "GREENS"
                    [Pref, XYref] = caching.cacheGreens(1,NeX,NeY,XY,...
                        path_ref_dir,acoustic_params,...
                        xy0_rel,bbox,boundary_type,Nw_greens);
                otherwise
                    error('Reference not supported')
            end
    
            assert(all(size(XYref) == size(XY), 'all'))

            fd = validation.roiFunctions(roi_type,bthickness,slice_source,bbox,xy0,r0);        
            Proi = validation.extractROI(P(:),XY,fd);
            Proi_ref = validation.extractROI(Pref(:),XY,fd);
            errs(i) = norm(abs(Proi) - abs(Proi_ref),2);
        end
    end
    
    errs_all{k} = errs;
    cpu_timings_all{k} = cpu_timings/num_avg_computations;
end

path_write_plot_f = @(filetype) sprintf('%s/conv_CPU_%s%0.1fX%0.1f_ref%s_xy%0.1fX%0.1f_f%0.1fhz_all.%s',...
            plot_path,sim_method_tuple{1},lx,ly,ref_method,xy0(1),xy0(2),f,filetype);

fig1 = figure(2);
plotting.plotConvergences(errs_all, cpu_timings_all, Porders_plot, 'CPU time [s]', legend_info, ...
    sprintf('Convergence\n%0.1fm x %0.1fm, f = %0.1fHz, ', lx, ly, f))

saveas(fig1,path_write_plot_f('png'))
saveas(fig1,path_write_plot_f('fig'))
    
function [xinputs,xlabelstr,plot_title,Porder] = setupSimulationParameters(sim_method_tuple)
    sim_method = sim_method_tuple{1};
    switch sim_method
        case "GREENS"            
            xinputs = 4:100:1004;
            xlabelstr = "Number wave functions (in total)";  
            plot_title = sprintf('%s convergence', sim_methods);
            Porder = {};
        case "WBM"            
            xinputs = 15:5:151;
            xlabelstr = "Number wave functions (in total)";
            plot_title = sprintf('%s convergence', sim_method);
            Porder = {};
        case "SEM"
            switch sim_method_tuple{2}
                case 1
                    Porder = 1;
                    xinputs = [20 40 60 80];        
                case 2
                    Porder = 2;
                    xinputs = [8 16 24 32];
                case 4
                    Porder = 4;
                    xinputs = [6 10 18 24];
                case 6
                    Porder = 6;
                    xinputs = [2 4 8 16];
                case 8
                    Porder = 8;
                    xinputs = [2 6 8 12];
            end            
            xlabelstr = "Number of elements";
            plot_title = sprintf('%s P=%i convergence', sim_method, Porder);
    end
end

function [Porder] = findMaxPolynomialOrder(NeX,NeY,max_num_nodes)
    Porder = 1;
    for i=1:50
        if (NeX*i+1)*(NeY*i+1) > max_num_nodes
            break;
        end
        Porder = i;
    end
end

% boundary source
% line_vel = [0 0; 0 0; 1 0; 0 0];
% line_press  = [1.0, 0.0, 0, 0];
% bound_velocities = [0 0; 1 0; 0 0; 0 0];
% bound_pressures = [0.0, 1.0, 0.0, 0.0];
% source_line_vel = models.SourceModel(models.SourceType.Line, line_vel);
% source_line_pres = models.SourceModel(models.SourceType.Line, line_press);
% domain_line_vel  = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, models.BoundaryCondition.Velocity, source_line_vel);
% domain_line_pres = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, models.BoundaryCondition.Pressure, source_line_pres);