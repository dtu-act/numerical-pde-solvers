function [p1,p2,sim1,sim2,iter] = run_FOURIER_FOURIER(tmax, xminmax, x0_pos, dx1, dx2, ...
    fn_source, source_partition, c, rho, scheme_order)
    import models.types.*
    import models.*
    import utilsDD.*
    
    cfl = sqrt(3);    
    
    l = xminmax(2) - xminmax(1);
    nmodes1 = l/(2*dx1);
    nmodes2 = l/(2*dx2);
    
    if ~(mod(dx1,dx2) < eps('double') || mod(dx2,dx1) < eps('double'))
        disp('WARNING: dx1 and dx2 are not a multiple of each other. Time resolution is set to the same in both partitions.')
        t_refine = 1;
    elseif abs(dx1 - dx2) < eps('double')
        t_refine = 1;
    else
        t_refine = 2;
    end
    
    dt2 = min(dx1,dx2)/(c*cfl);
    dt1 = dt2*t_refine;    
    
    if dt1 > dx1/(c*cfl) || dt2 > dx2/(c*cfl)
        error('dt not satisfying CFL!')
    end
    
    iter = tmax/dt1;
    if rem(iter,1) > eps('single')
        iter = ceil(iter);
        dt1_adjust = tmax/iter;
        fprintf('WARNING: dt1 adjusted from %e to %e to match tmax\n', dt1, dt1_adjust)
        dt1 = dt1_adjust;
        dt2 = dt1/t_refine;
    end
    
    x1d_1 = xminmax(1):dx1:xminmax(2);
    x1d_2 = xminmax(1):dx2:xminmax(2);
    
    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax, BoundaryType.Neumann);
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax, BoundaryType.Neumann);
    
    sourceFactor = sourceScalingFourier()*2;
    if source_partition == "LEFT"
        src1 = defaultSource(fn_source,sourceFactor,domain1,x0_pos);
        src2 = {};
    elseif source_partition == "RIGHT"
        src2 = defaultSource(fn_source,sourceFactor,domain2,x0_pos);
        src1 = {};
    else
        error('No such domain')
    end
        
    sim1 = Simulation1D(SolverType.FOURIER, domain1, src1, solvers.CustomFourier(nmodes1));
    sim2 = Simulation1D(SolverType.FOURIER, domain2, src2, solvers.CustomFourier(nmodes2));
    
    [p1,p2] = solvers.runSolver(iter, sim1, sim2, scheme_order);
end