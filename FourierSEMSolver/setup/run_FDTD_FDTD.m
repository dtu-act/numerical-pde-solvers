function [p1,p2,sim1,sim2,iter] = run_FDTD_FDTD(tmax, xminmax, dx1, dx2, fn_source, c, rho, scheme_order)    
    cfl = sqrt(3);
    
    x1d_1 = xminmax(1):dx1:xminmax(2);
    x1d_2 = xminmax(1):dx2:xminmax(2);
    
    dt1 = dx1/(c*cfl);
    dt2 = dt1/2; %dx2/(c*cfl);
    
    iter = ceil(tmax/min(dt1,dt2));
    
    sourceFactor = sourceScalingFDTD(scheme_order)*8;
    
    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1);
    sim1 = Simulation1D(SolverType.FDTD, domain1, src1, CustomFDTD(scheme_order));
        
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax, BoundaryType.Neumann);    
    sim2 = Simulation1D(SolverType.FDTD, domain2, {}, CustomFDTD(scheme_order));
    
    [p1,p2] = runSolver(iter, sim1, sim2, scheme_order);
end