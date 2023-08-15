function [state] = solverWE1D(state)
    import models.types.*
    import solvers.*

    switch state.solver_type
        case SolverType.FDTD
            solverWE1D_FDTD(state);
        case SolverType.FOURIER
            solverWE1D_FOURIER(state);
        case SolverType.SEM
            solverWE1D_SEM(state);
        otherwise
            error('Solver not supported')
    end    
end