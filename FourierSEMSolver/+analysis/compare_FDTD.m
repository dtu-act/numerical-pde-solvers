function compare_FDTD(iter, xminmax, dx, fn_source, c, rho, scheme_order, interface_order)    
    dt = dx/(c*sqrt(3));
    sourceFactor = sourceScalingFDTD(scheme_order);
    
    x1d = xminmax(1):dx:xminmax(2);
    
    domain1 = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1);
    state1 = Simulation1D(SolverType.FDTD, domain1, src1, CustomFDTD(scheme_order));
    
    run(iter, state1, interface_order, true);
end