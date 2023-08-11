function compare_SEM_SEM(adjust, iter, xminmax, dx, fn_source, c, rho, P1_order, ...\
        P2_order, interface_order1, interface_order2, iface1, iface2)
    
    x0 = (xminmax(2)-xminmax(1))/2;

    xminmax_1 = xminmax;
    xminmax_2 = xminmax;
    
    if iface1.Nk > 0        
        xminmax_1(2) = xminmax_1(2)+iface1.dx*2*iface1.Nk;
    end
    
    if iface2.Nk > 0       
        xminmax_2(2) = xminmax_2(2)+iface2.dx*2*iface2.Nk;
    end
    
    cfl = sqrt(3);
    [x1d1,custom1,dxMin1] = setupSEMCustom(xminmax_1, dx, P1_order, iface1);
    [x1d2,custom2,dxMin2] = setupSEMCustom(xminmax_2, dx, P2_order, iface2);
    
    dt1 = min(dxMin1,dxMin2)/(c*cfl);
    dt2 = dt1; %dxMin2/(c*cfl);
    
    sourceFactor1 = sourceScalingSEM(1)*8;
    sourceFactor2 = sourceScalingSEM(1)*8;
    
    domain1 = Domain1D(x1d1, dx, dt1, c, rho, xminmax_1, BoundaryType.Neumann);
    src1 = defaultSource(fn_source,sourceFactor1,domain1,x0);
    sim1 = Simulation1D(SolverType.SEM, domain1, src1, custom1);
    
    domain2 = Domain1D(x1d2, dx, dt2, c, rho, xminmax_2, BoundaryType.Neumann);
    src2 = defaultSource(fn_source,sourceFactor2,domain2,x0);
    sim2 = Simulation1D(SolverType.SEM, domain2, src2, custom2);
    
    run_compare(adjust, iter, sim1, sim2, "none", interface_order1, interface_order2);
end