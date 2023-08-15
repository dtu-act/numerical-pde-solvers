function [p1,p2,sim1,sim2,iter] = run_SEM_SEM(tmax, xminmax, x0_pos, dx_1, dx_2, fn_source, c, rho, P_order, interface_order, sem_iface)
    cfl = sqrt(3);
    
    if sem_iface.Nk > 0
        L_iface = sem_iface.dx*2*sem_iface.Nk;
        xminmax(2) = xminmax(2)+L_iface;
    end    

    [x1d_1,custom_1,dxMin_1] = setup.setupSEMCustom(xminmax, dx_1, P_order, sem_iface);
    [x1d_2,custom_2,dxMin_2] = setup.setupSEMCustom(xminmax, dx_2, P_order, sem_iface);
    
    dt_1 = dxMin_1/(c*cfl);
    dt_2 = dxMin_2/(c*cfl);
    
    iter = tmax/dt_1;
    
    sourceFactor = sourceScalingSEM(P_order);
    
    domain1 = Domain1D(x1d_1, dx_1, dt_1, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1,x0_pos);
    sim1 = Simulation1D(SolverType.SEM, domain1, src1, custom_1);
        
    domain2 = Domain1D(x1d_2, dx_2, dt_2, c, rho, xminmax, BoundaryType.Neumann);
    sim2 = Simulation1D(SolverType.SEM, domain2, {}, custom_2);
    
    [p1,p2] = runSolver(iter, sim1, sim2, interface_order);
end