function [P, XY, conn] = readWriteSEM(path,Porder,NeX,NeY,bbox,boundary_type,xy0,aparams)
    if isfile(path)
        fprintf('[cache] Loading SEM solution from disk ... \n')
        data = load(path,'P','XY','conn');
        P = data.P;
        XY = data.XY;
        conn = data.conn;
    else
        fprintf('[cache] Calculating SEM solution ... \n')
        
        source = models.PointSource(xy0,1);
        
        solverSEM = solvers.HelmholtzSolver2D(aparams);
        solverSEM.setupMeshUniform(bbox,NeX,NeY);
        solverSEM.setupSolver(Porder, boundary_type, source);                
        [P, XY] = solverSEM.solve();
        conn = solverSEM.conn;
        
        fprintf('[cache] Writing SEM solution to disk: %s... \n', path)            
        save(path, 'P', 'XY', 'conn');
    end
end