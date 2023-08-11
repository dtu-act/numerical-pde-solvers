function [p1,p2,sim1,sim2,iter] = run_FOURIER_FDTD(tmax, xminmax, dx1, dx2, fn_source, c, rho, scheme_order)
    cfl1 = sqrt(3);
    cfl2 = sqrt(3);
    
    l = xminmax(2) - xminmax(1);
    nmodes1 = l/(2*dx1);
    
    dt1 = dx1/(c*cfl1);
    dt2 = dt1/2; %dx2/(c*cfl2);
    
    iter = ceil(tmax/min(dt1,dt2));
    
    x1d_1 = xminmax(1):dx1:xminmax(2);
    x1d_2 = xminmax(1):dx2:xminmax(2);
    
    sourceFactor = sourceScalingFourier();
    
    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1);
    sim1 = Simulation1D(SolverType.FOURIER, domain1, src1, CustomFourier(nmodes1));
        
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax, BoundaryType.Neumann);
    sim2 = Simulation1D(SolverType.FDTD, domain2, {}, CustomFDTD(scheme_order));
    
    [p1,p2] = runSolver(iter, sim1, sim2, scheme_order);
end