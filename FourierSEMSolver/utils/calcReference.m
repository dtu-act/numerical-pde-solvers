function [p, s1] = calcReference(xminmax, ppw, dt, fmax, x0_pos, c, rho, iter, tnorm)    
    [dx,nmodes] = calcSpatial(xminmax(2),fmax,ppw,c);
    if dt > dx/c
        error('dt not satisfying CFL!')
    end
    
    x1d = xminmax(1):dx:xminmax(end);
    assert(abs(x1d(end) - xminmax(end)) < eps('double'))
    
    sourceFactor = sourceScalingFourier()*2;
    
    domain1 = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
    src1 = defaultSource(fmax,sourceFactor,domain1,x0_pos);
    assert(abs(src1.x0 - x0_pos) < eps('single'))
    
    s1 = Simulation1D(SolverType.FOURIER, domain1, src1, CustomFourier(nmodes));
    
    p = runSingleDomainSolver(iter, s1);
    
    n_norm = round(tnorm/dt);
    norm_factor = max(p(n_norm,:));
    
    p = p/norm_factor;
end