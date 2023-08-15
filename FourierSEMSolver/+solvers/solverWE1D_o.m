function [p,n] = solverWE1D_o(n,solver_type,p,xminmax,c,dt,custom,F,bound_type)
    import utilsDD.*
    import models.types.*
    import solvers.*

    switch solver_type
        case SolverType.FDTD
            error('NOT IMPLEMENTED (use non-optimized version instead)')
            [pc,pp] = solverWE1D_FDTD(state);
        case SolverType.FOURIER
            [p,n] = solverWE1D_FOURIER_o(n,p,xminmax,c,dt,custom.nmodes,F,bound_type);            
        case SolverType.SEM
            [p,n] = solverWE1D_SEM_o(n,p,c,dt,custom.Mx,custom.Sx,F,bound_type);
        otherwise
            error('Solver not supported')
    end    
end