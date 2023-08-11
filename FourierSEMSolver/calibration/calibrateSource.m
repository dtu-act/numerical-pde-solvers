function [sourceFactor] = calibrateSource(method, dx, nmodes, c, rho, xminmax, fn_source, order, iter)
    sourceFactor = 1;
    
    cfl = sqrt(3);    
    
    switch method
        case SolverType.FDTD
            method_str = 'FDTD';            
            x1d = xminmax(1):dx:xminmax(2);
            dt = dx/(c*cfl); % use CFL from FDTD schemes
            domain = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
            src = defaultSource(fn_source,sourceFactor,domain);        
            sim = Simulation1D(SolverType.FDTD, domain, src, CustomFDTD(order));
        case SolverType.FOURIER            
            method_str = 'FOURIER';
            x1d = xminmax(1):dx:xminmax(2);
            dt = dx/(c*cfl); % use CFL from FDTD schemes
            domain = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
            src = defaultSource(fn_source,sourceFactor,domain);        
            sim = Simulation1D(SolverType.FOURIER, domain, src, CustomFourier(nmodes));
        case SolverType.SEM
            %if order == 1
                disp('Time stepping refined by 2')
                ts_factor = 2;
            %else
            %    ts_factor = 1;
            %end
            method_str = 'SEM';
            [x1d,custom,dxMin] = setupSEMCustom(xminmax, dx, order);
            dt = dxMin/(c*cfl)/ts_factor;
            domain = Domain1D(x1d, dx, dt, c, rho, xminmax, BoundaryType.Neumann);
            src = defaultSource(fn_source,sourceFactor,domain);
            sim = Simulation1D(SolverType.SEM, domain, src, custom);
        otherwise
            error('method not supported')
    end    
    
    while sim.n < iter              
        sim = solverWE1D(sim);
        
        t = sim.n*sim.domain.dt;
        
        figure(1)
        plot(domain.x1d, sim.p_current,'-r')    
        xlabel('n')
        ylabel('p')
        title(sprintf('%s time = %.4f [s]', method_str, t))
        
        sourceFactor = 0.5/max(sim.p_current);
    end
end