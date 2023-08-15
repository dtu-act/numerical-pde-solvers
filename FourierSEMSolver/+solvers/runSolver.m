function [p1_all,p2_all] = runSolver(iter, sim1, sim2, interface_order)
    import solvers.*
    import interpolation.*

    doTemporalAdjustment = abs(sim1.domain.dt-sim2.domain.dt) > eps('double');
    
    p1_all = zeros(iter,length(sim1.p_current));
    p2_all = zeros(iter,length(sim2.p_current));
    
    %figure(3);
    %n_stop = 100;
    
    tstart_tot = tic;
    
    tot_time1 = 0;
    tot_time2 = 0;
    tot_time3 = 0;
    tot_time4 = 0;
    tot_time5 = 0;
    
    for n=1:iter
        if doTemporalAdjustment
            if sim1.domain.dt < sim2.domain.dt
                sim1_tmp = sim1;                
                sim1 = sim2;
                sim2 = sim1_tmp;
            end
            
            tstart = tic;
            %% n-1 |+| n -> n+1
            % SOLVE (FULL) STATE 1
            solverWE1D(sim1);

            %% n-1/2 |+| n -> n+1/2
            % SOLVE STATE 2
            state2_1 = solverWE1D(sim2);
            tot_time1 = tot_time1 + toc(tstart);
            
%             if n > n_stop
%                 x1d_1 = sim1.domain.x1d;
%                 x1d_2 = sim1.domain.x1d(end) - state2_1.domain.x1d(1) + state2_1.domain.x1d;
% 
%                 plot(x1d_1,sim1.p_current,'-o')
%                 hold on
%                 plot(x1d_2,state2_1.p_current,'-o')
%                 hold off
%                 ylim([-1,1])
%                 drawnow
%             end
            
            tstart = tic;
            % INTERFACE handling, compensating using (ONLY) PREVIOUS time step n
            % sim 1: n+1   (prev is n)
            % sim 2: n+1/2 (prev is n)
            % NOTE: compensate for interface pressures BEFORE temporal interpolation           
            [sim1, state2_1_i] = solverWE1DInterface(sim1, state2_1, interface_order);
            tot_time2 = tot_time2 + toc(tstart);
            
            %% n |+| n+1/2 -> n+1
            state2_2 = solverWE1D(state2_1_i);

            tstart = tic;
            % create prev n+1/2 for sim 1 by interpolation
            [state1_1_i] = adjustTemporal(sim1,sim2.domain.dt);
            tot_time3 = tot_time3 + toc(tstart);
            
            tstart = tic;
            % adjust previous and update dt
            state1_2_i = copy(sim1);
            state1_2_i.domain.dt = sim2.domain.dt;
            state1_2_i.adjust_prev(state1_1_i.p_prev);
            tot_time4 = tot_time4 + toc(tstart);
            
            tstart = tic;
            % INTERFACE handling, compensating using (ONLY) PREVIOUS time step n+1/2
            % sim 1: n+1 (prev is n+1/2)
            % sim 2: n+1 (prev is n+1/2)
            [~, sim2] = solverWE1DInterface(state1_2_i, state2_2, interface_order);
            tot_time5 = tot_time5 + toc(tstart);
        else
            %% n-1 |+| n -> n+1
            % SOLVE (FULL) STATE 1
            solverWE1D(sim1);
        
            %% n-1 |+| n -> n+1
            % SOLVE STATE 2
            solverWE1D(sim2);
            
            % INTERFACE handling            
            [sim1, sim2] = solverWE1DInterface(sim1, sim2, interface_order);
        end
        
        p1_all(n,:) = sim1.p_current;
        p2_all(n,:) = sim2.p_current;
        
        if mod(n,100) == 0
            fprintf('Calculating n=%i ...\n',n)
        end
    end
    
    tot_time = toc(tstart_tot);
    
    fprintf('TOTAL: %0.3f \n', tot_time)
    fprintf('solverWE1D: %0.3f \n', tot_time1)
    fprintf('solverWE1DInterface: %0.3f \n', tot_time2)
    fprintf('adjustTemporal: %0.3f \n', tot_time3)
    fprintf('solverWE1DInterface: %0.3f \n', tot_time5)  
end