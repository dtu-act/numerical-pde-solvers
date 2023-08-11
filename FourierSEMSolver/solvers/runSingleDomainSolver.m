function [p1_all] = runSingleDomainSolver(iter, sim)
    p1_all = zeros(iter,length(sim.p_current));
    
    for n=1:iter
        % SOLVE (FULL) STATE 1
        sim = solverWE1D(sim);
        
        p1_all(n,:) = sim.p_current;
        
%         if mod(n,100) == 0
%             fprintf('Calculating n=%i ...\n',n)
%         end
    end
end