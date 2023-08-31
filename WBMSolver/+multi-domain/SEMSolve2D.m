function [PSEM, XY] = SEMSolve2D(P_order, NeX, NeY, bbox, xy0, boundary_type, source_type, acoustic_params)
    mesh_type = "uniform";
    h0 = 0.05;

    % FORCING
    switch source_type
        case models.SourceType.PointSource
            source = models.PointSource(xy0,-1);
        case models.sourceType.Gaussian
            sigma = 0.2; % -> 1000Hz
            source = setup.gauss2DSetup(c,sigma,xy0);
        case models.sourceType.Gaussian
            source = models.SourceModelNone();
        otherwise
            error('source type not supported')
    end

    % MESH
    bound_detect_f = meshing.createBoundRectDetectFunc(bbox, 1e-10);
    bound_cond_f = @(x,y) 0; % only used for pressure boundary conditions

    switch mesh_type
        case 'uniform'
            mesh = meshing.mesh2D(bbox,NeX,NeY);
            etov = meshing.createEToV(NeX,NeY);
            domain = sem.Domain2D(bbox,mesh,source,etov,boundary_type,bound_detect_f,bound_cond_f);
        case 'nonuniform'        
            [mesh, etov] = meshing.distMeshGeneration(bbox,xy0,h0);
            domain = sem.Domain2D(bbox,mesh,source,etov,boundary_type,bound_detect_f,bound_cond_f);
        otherwise
            error('mesh type not supported')
    end

    simulation_params = sem.SimulationParameters2D(P_order);

    %==== CALCULATE =====%
    [PSEM, XY] = solvers.solveHelmholtz2D(domain,acoustic_params.k,simulation_params);
end