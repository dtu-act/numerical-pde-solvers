function compare_SEM(iter, xminmax, dx, fn_source, c, rho, P_order, interface_order)
    sourceFactor = sourceScalingSEM(P_order)/2;
    
    [x1d,custom1,dxMin] = setupSEMCustom(xminmax, dx, P_order);
    
    cfl = sqrt(3);
    dt = dxMin/(c*cfl)/2;
    
    domain1 = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1);
    sim1 = Simulation1D(SolverType.SEM, domain1, src1, custom1);
    
    run(iter, sim1, interface_order, true);
end