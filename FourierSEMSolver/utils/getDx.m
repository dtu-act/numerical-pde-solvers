function dx = getDx(sim)
    if sim.solver_type == SolverType.SEM
        if  isempty(sim.custom.interface)
            if sim.custom.scheme_order ~= 1
                error('NOT IMPLEMENTED')
            end
            dx = sim.domain.dx;
        else
            dx = sim.custom.interface.dx;
        end
    else        
        dx = sim.domain.dx;
    end
end