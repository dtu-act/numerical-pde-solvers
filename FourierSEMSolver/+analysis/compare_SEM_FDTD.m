function compare_SEM_FDTD(iter, xminmax, dx, fn_source, c, rho, scheme_order1, scheme_order2, interface_order, iface)
    import models.types.*
    import models.*
    import setup.*
    import utilsDD.*

    x0 = (xminmax(2)-xminmax(1))/2;

    xminmax_1 = xminmax;
    
    if isempty(iface)
        iface = SEMUniformInterface(dx, 0, 1, false, InterfaceLocation1D.LEFT);
    end
    
    if iface.loc == InterfaceLocation1D.LEFT
        xminmax_1(1) = xminmax_1(1)-iface.dx*iface.P_order*iface.Nk;
    else
        xminmax_1(2) = xminmax_1(2)+iface.dx*iface.P_order*iface.Nk;
    end
    
    xminmax_2 = xminmax_1;
    
    cfl = sqrt(3);
    [x1d_1,custom1,dxMin] = setupSEMCustom(xminmax_1, dx, scheme_order1, iface);  
    
    dx1 = dx;
    dx2 = iface.dx;
    dt1 = 10^-5; %dxMin/(c*cfl);
    dt2 = dt1; %dx/(c*cfl);
    
    assert(dt1 < dxMin/(c*cfl))
    
    sourceFactor1 = sourceScalingSEM(1)/1;
    sourceFactor2 = sourceScalingFDTD(scheme_order2)/4;
    
    x1d_2 = xminmax_2(1):dx2:xminmax_2(2);
    
    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax_1, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor1,domain1,x0);
    sim1 = Simulation1D(SolverType.SEM, domain1, src1, custom1);
    
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax_2, BoundaryType.Neumann);
    src2 = defaultSource(fn_source,sourceFactor2,domain2,x0);
    sim2 = Simulation1D(SolverType.FDTD, domain2, src2, solvers.CustomFDTD(scheme_order2));
    
    analysis.run_compare(iter, sim1, sim2, "none", interface_order, interface_order, iface.loc);
end