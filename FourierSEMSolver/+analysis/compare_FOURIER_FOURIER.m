function compare_FOURIER_FOURIER(adjust, iter, xminmax, dx, nmodes, fn_source, c, rho, scheme_order1, scheme_order2, iloc)
    import utilsDD.*
    import models.types.*
    import models.*
    
    cfl = sqrt(3);
    x0 = (xminmax(2)-xminmax(1))/2;

    [dx1,dx2,dt1,dt2,nmodes1,nmodes2] = calcResolutionParameters(adjust, dx, c, cfl, nmodes);

    sourceFactor = sourceScalingFourier();

    x1d_1 = xminmax(1):dx1:xminmax(2);
    x1d_2 = xminmax(1):dx2:xminmax(2);

    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor,domain1,x0);
    sim1 = Simulation1D(SolverType.FOURIER, domain1, src1, solvers.CustomFourier(nmodes1));
    
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax, BoundaryType.Neumann);
    src2 = defaultSource(fn_source,sourceFactor,domain2,x0);
    sim2 = Simulation1D(SolverType.FOURIER, domain2, src2, solvers.CustomFourier(nmodes2));
    
    analysis.run_compare(iter, sim1, sim2, adjust, scheme_order1, scheme_order2, iloc)
end