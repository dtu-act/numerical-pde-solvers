function compare_FDTD_FDTD(adjust, iter, xminmax, dx, fn_source, c, rho, scheme_order1, scheme_order2, iloc)
    cfl = sqrt(3);
    x0 = (xminmax(2)-xminmax(1))/2;

    [dx1,dx2,dt1,dt2] = calcResolutionParameters(adjust, dx, c, cfl);
    
    sourceFactor1 = sourceScalingFDTD(scheme_order1);
    if dx1 ~= dx2
        sourceFactor2 = sourceScalingFDTD(scheme_order2)/4;
    else
        sourceFactor2 = sourceScalingFDTD(scheme_order2);
    end
    
    x1d_1 = xminmax(1):dx1:xminmax(2);
    x1d_2 = xminmax(1):dx2:xminmax(2);
    
    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor1,domain1,x0);
    sim1 = Simulation1D(SolverType.FDTD, domain1, src1, CustomFDTD(scheme_order1));
    
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax, BoundaryType.Neumann);
    src2 = defaultSource(fn_source,sourceFactor2,domain2,x0);
    sim2 = Simulation1D(SolverType.FDTD, domain2, src2, CustomFDTD(scheme_order2));
    
    run_compare(iter, sim1, sim2, adjust, scheme_order1, scheme_order2, iloc);
end