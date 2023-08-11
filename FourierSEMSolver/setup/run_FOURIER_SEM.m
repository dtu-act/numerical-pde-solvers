function [p1,p2,sim1,sim2,iter] = run_FOURIER_SEM(tmax, xminmax_1, xminmax_2, ...\
    dx1, dx2, fmax, source_partition, x0_pos, c, rho, P_order, interface_order, iface)
    
    run_optimized = true;

    cfl = sqrt(3);
    
    l = xminmax_1(2) - xminmax_1(1);
    nmodes1 = l/(2*dx1);
    
    [x1d_2,custom2,dxMin2] = setupSEMCustom(xminmax_2, dx2, P_order, iface);
    
    if ~(mod(dx1,dx2) < eps('double') || mod(dx2,dx1) < eps('double'))
        disp('WARNING: dx1 and dx2 are not a multiple of each other. Time resolution is set to the same in both partitions.')
        t_refine = 1;
    elseif abs(dx1 - dx2) < eps('double')
        t_refine = 1;
    else
        t_refine = 2;
    end
    
    dt2 = min(dx1,dxMin2)/(c*cfl);
    dt1 = dt2*t_refine;
    
    iter = tmax/dt1;
    if rem(iter,1) > eps('single')
        iter = ceil(iter);
        dt1_adjust = tmax/iter;
        fprintf('WARNING: dt1 adjusted from %e to %e to match tmax\n', dt1, dt1_adjust)
        dt1 = dt1_adjust;
        dt2 = dt1/t_refine;
    end
    
    if dt1 > dx1/(c*cfl) || dt2 > min(dx1,dxMin2)/(c*cfl)
        error('dt not satisfying CFL!')
    end
        
    x1d_1 = xminmax_1(1):dx1:xminmax_1(2); 
    
    domain1 = Domain1D(x1d_1, dx1, dt1, c, rho, xminmax_1, BoundaryType.Neumann);    
    domain2 = Domain1D(x1d_2, dx2, dt2, c, rho, xminmax_2, BoundaryType.Neumann);    
    
    if source_partition == "LEFT"
        sourceFactor = sourceScalingFourier()*2;
        src1 = defaultSource(fmax,sourceFactor,domain1,x0_pos);
        src2 = {};
    elseif source_partition == "RIGHT"
        sourceFactor = sourceScalingSEM(P_order)*8;
        src2 = defaultSource(fmax,sourceFactor,domain2,x0_pos);
        src1 = {};
    else
        error('No such domain')
    end
    
    sim1 = Simulation1D(SolverType.FOURIER, domain1, src1, CustomFourier(nmodes1));
    sim2 = Simulation1D(SolverType.SEM, domain2, src2, custom2);
    
    if run_optimized
        [p1,p2] = runSolver_o(iter, sim1, sim2, interface_order);
    else
        [p1,p2] = runSolver(iter, sim1, sim2, interface_order);
    end
end