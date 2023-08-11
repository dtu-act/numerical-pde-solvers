function [state] = FDTDSimulationTestSetup(x1, domain, source, custom, nhistory)
    state = Simulation1D(SolverType.FDTD, domain, source, custom);
    N = length(state.p_current);
    
    p_history = zeros(N,nhistory);
    alphas = linspace(1,0.5,nhistory);
    
    for i=1:nhistory
        p_history(:,nhistory-(i-1)) = alphas(i)*cos(x1);
    end
    
    state = state.setHistory(p_history);
end