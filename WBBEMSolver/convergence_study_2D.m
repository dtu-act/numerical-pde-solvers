clear all
close all
addpath('../shared')
addpath('../SEMSolvers')

[path_ref_dir, plot_path] = paths.setupPaths("LOCAL");

%% SIMULATION METHOD/REFERENCE
sim_method = "WBM"; % "GREENS", "WBM", "SEM"
ref_method = "SEM"; % "GREENS" | "SEM"

Porder_ref = 6;   % used only when comparing with SEM
Nw_greens = 4000; % used only when comparing with Greens

% SETUP
%roi_types = ["fd_full", "fd_boundaries", "fd_axial", "fd_bsourcenoaxial", "fd_fullnobounds", "fd_Q1Q2Q3Q4"];
%legend_info = {'Full', 'Boundaries', 'Axial', 'Circle around source', 'Inner w/o boundaries', 'Quadrants w/o axial/boundaries'};
%line_styles = ["-ok", "-+g", "-*r", "-xm", "-dc", "-^b"];

roi_types = ["fd_full", "fd_fullnobounds"];
legend_info = {'Full domain', 'Inner domain'};
line_styles = ["-ob", "-or"];

bthickness = 0.01;
slice_source = 0.2;

% acoustic parameters
f = 300;
c = 343.0;
rho = 1.225;
acoustic_params = models.AcousticParameters(f,c,rho);
boundary_type = models.BoundaryCondition.Velocity;

% domain size
lx = 2;
ly = 2;
bbox = [0,0;lx,ly];
xy0_rel = [3/5, 3/5]; % source relative location

grid_res_ref = 0.1;
r0 = 0.4; % distance from source to neglect
    
% setup resolutions and parameters for convergence test
switch sim_method
    case "GREENS"
        Porders_theoretical_plot = [0];                
        sim_resolutions = 4:100:1004;
        sim_label = "Number of wave functions (total)";  
        convplottitle = sprintf('%s convergence', sim_method);
    case "WBM"
        Porders_theoretical_plot = [1,3];
        sim_resolutions = 15:5:151;
        sim_label = "Number of wave functions (total)";
        convplottitle = sprintf('%s convergence', sim_method);
    case "SEM"
        Porders_theoretical_plot = [2,4];
        solverSEM = solvers.HelmholtzSolver2D(acoustic_params);
        
        Porder = 4;
        sim_resolutions = [6 10 18 24];
        
        sim_label = "Number of elements";
        convplottitle = sprintf('%s P=%i convergence', sim_method, Porder);
    otherwise
        error('simulation method not supported')
end

disp('*****Simulation info*****')
fprintf('Domain size: %0.1fx%0.1f\n', lx,ly)
fprintf('Nw Greens: %i\n', Nw_greens)
fprintf('s_{xy} = (%0.2f,%0.2f)\n\n', lx*xy0_rel(1), ly*xy0_rel(2))

% containers
errs_all = cell(length(roi_types),1);
errs = zeros(1,length(sim_resolutions));
sim_resolutions_out = zeros(1,length(sim_resolutions));

%   START SIMULATING    %

for k=1:length(roi_types)
    roi_type = roi_types(k);
    
    for i=1:length(sim_resolutions)
        xinput = sim_resolutions(i);
        
        switch sim_method            
            case "GREENS"
                [xy0, ~, XY] = setupDomain(lx,ly,xy0_rel,boundary_type,grid_res_ref);
                
                NeX = lx/grid_res_ref; NeY = ly/grid_res_ref;
                
                fprintf('*** Compute Greens for Nw=%i...\n', xinput)
                assert(lx == ly)
                P = greens.greensNeumannUnstructuredSpatial2D(xinput,xinput,XY,xy0,acoustic_params.k);
                sim_resolutions_out(i) = xinput*xinput;
            case "WBM"
                [xy0, domain] = setupDomain(lx,ly,xy0_rel,boundary_type,grid_res_ref);
                
                NeX = lx/grid_res_ref; NeY = ly/grid_res_ref;
                
                fprintf('*** Compute WBM for w=%i...\n', xinput)
                [P, XY, wavecount] = wbm.solveHelmholtzSingle2D(domain, acoustic_params, xinput, xinput, grid_res_ref);
                sim_resolutions_out(i) = wavecount.nw;
            case "SEM"
                NeX = xinput; NeY = xinput;
                
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
        
                sim_resolutions_out(i) = xinput;
            otherwise
                error('Simulation method not supported')
        end
        
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
    errs_all{k} = errs;
    
    % PLOT     
    fig = figure(2);

    plotting.plotConvergenceMesh(errs,sim_resolutions_out,Proi,XY,xy0,sim_label,sim_method, ...
        convplottitle,"-b",45);

    roi_type_str = split(roi_type, '_');
    roi_type_str = roi_type_str(2);
    
    path_write_plot = sprintf('%s/conv_%s%0.1fX%0.1f_ref%s_xy%0.1fX%0.1f_f%0.1fhz_%s.png',...
            plot_path,sim_method,lx,ly,ref_method,xy0(1),xy0(2),f,roi_type_str);
    
    titlestr = sprintf('%s, %im x %im, s_{xy}=(%0.1f, %0.1f)\nROI: %s', ...
        convplottitle,lx,ly,xy0(1),xy0(2), roi_type_str);
    
    sgtitle(titlestr)
    saveas(fig,path_write_plot)    
    close(fig)
end

path_write_plot = sprintf('%s/conv_%s%0.1fX%0.1f_ref%s_xy%0.1fX%0.1f_f%0.1fhz_all.png',...
            plot_path,sim_method,lx,ly,ref_method,xy0(1),xy0(2),f);

fig = figure(3);
xvalues_all = cell(size(errs_all,1),1);
[xvalues_all{:}] = deal(sim_resolutions_out);
title_str = ''; %sprintf('%s convergences for ROIs\n%0.1fm x %0.1fm, f = %0.1fHz, ', sim_method, lx, ly, f);

plotting.plotConvergences(errs_all, xvalues_all, Porders_theoretical_plot, sim_label, legend_info, ...,
    line_styles, title_str)

saveas(fig,path_write_plot)

function [xy0, domain, XY] = setupDomain(lx,ly,xy0_rel,boundary_type,grid_res)
    xy0 = [lx,ly].*xy0_rel;
        
    gcoord = [0.0 0.0; lx*1.0 0.0; lx*1.0 ly*1.0; 0.0 ly*1.0]; % coordinates at boundaries
    origin = [0,0];
    bound_nodes = [1 2; 2 3; 3 4; 4 1];

    source = wbm.PointSource(xy0,1);
    domain = wbm.Domain2D(lx, ly, gcoord, bound_nodes, origin, boundary_type, source);
    
    [X,Y] = meshgrid(0:grid_res:lx, 0:grid_res:ly);
    xs = domain.origin(1) + X; % include offset -> original position
    ys = domain.origin(2) + Y; % include offset -> original position
    XY = [xs(:),ys(:)];
end