function [p1,p2,sim1,sim2,iter] = run_SEM_SEM(tmax, xminmax, dx, fn_source, c, rho, P_order, interface_order, sem_iface)    
    cfl = sqrt(3);
    
    if sem_iface.Nk > 0
        L_iface = sem_iface.dx*2*sem_iface.Nk;
        xminmax(2) = xminmax(2)+L_iface;
    end
        
    [x1d,custom,dxMin] = setupSEMCustom(xminmax, dx, P_order, sem_iface);
    
    dt = dxMin/(c*cfl);
    
    iter = tmax/dt;
    
    sourceFactor = sourceScalingSEM(P_order);
    
    domain1 = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1);
    sim1 = Simulation1D(SolverType.SEM, domain1, src1, custom);
        
    domain2 = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
    sim2 = Simulation1D(SolverType.SEM, domain2, {}, custom);
    
    [p1,p2] = runSolver(iter, sim1, sim2, interface_order);
end