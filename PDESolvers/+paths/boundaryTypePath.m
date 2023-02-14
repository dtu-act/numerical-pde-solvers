function path_dir = boundaryTypePath(path_dir,boundary_type)
        
    switch boundary_type
        case models.BoundaryCondition.Pressure
            boundstr = 'dirichlet';
        case models.BoundaryCondition.Velocity
            boundstr = 'neumann';
        case models.BoundaryCondition.Impedance
            boundstr = 'impedance';
        otherwise
            error('boundary type not supported')
    end
    
    path_dir = sprintf('%s/%s',path_dir,boundstr);
end