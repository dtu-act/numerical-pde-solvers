function [p1_all,p2_all] = runSolver_o(iter, sim1, sim2, order)    
    doTemporalAdjustment = abs(sim1.domain.dt-sim2.domain.dt) > eps('double');
    
    p1_all = zeros(iter,length(sim1.p_current));
    p2_all = zeros(iter,length(sim2.p_current));
    
    % shared
    c = sim1.domain.c;    
    
    % left domain
    solver_type1 = sim1.solver_type;
    p1 = sim1.p_history;
    dt1 = sim1.domain.dt;
    dx1 = getDx(sim1);
    custom1 = sim1.custom;
    F1 = sim1.src.F;
    bound_type1 = sim1.domain.boundary;
    xminmax1 = [sim1.domain.xmin,sim1.domain.xmax];
    x1d1 = sim1.domain.x1d;
    
    % right domain
    solver_type2 = sim2.solver_type;
    p2 = sim2.p_history;
    dt2 = sim2.domain.dt;
    dx2 = getDx(sim2);
    custom2 = sim2.custom;
    F2 = sim2.src.F;
    bound_type2 = sim2.domain.boundary;
    xminmax2 = [sim2.domain.xmin,sim2.domain.xmax];
    x1d2 = sim2.domain.x1d;
    
    if doTemporalAdjustment
        % time interpolation for sim1
        p1_history = sim1.p_history;
        steps_in = size(p1_history,2); 
        steps_out = 2;
        x = 0:dt1:(steps_in-1)*dt1; % number of steps of timesteps (history)
        xx = 0:dt2:(steps_out-1)*dt2;
    end
    
    N_iface = order/2 + 1;
    x1d1_dx2 = fliplr(xminmax1(2):-dx2:(xminmax1(2) - dx2*(N_iface-1))); % flip to ensure cont. at right interface    
    x1d2_dx1 = xminmax2(1):dx1:(xminmax2(1) + dx1*(N_iface-1));
    %x1d1_dx2 = fliplr(xminmax1(2):-dx2:(xminmax1(1)));    
    %x1d2_dx1 = xminmax2(1):dx1:(xminmax2(2));; % flip to ensure cont. at right interface    
        
    tot_time_solver1 = 0;
    tot_time_solver2 = 0;
    tot_time_interface = 0;
    tot_time_tspline = 0;
    
    if dt1 < dt2
        error("NOT IMPLEMENTED: dt1 < dt2")
    end
    
    n1 = 0;
    n2 = 0;
    
    tstart_tot = tic;
    
    fprintf('SOLVER LOOP STARTING ...\n')
    
    for n=1:iter
        if doTemporalAdjustment                        
            tstart = tic;
            %% n-1 |+| n -> n+1
            % SOLVE (FULL) STATE 1
            [p1,n1] = solverWE1D_o(n1,solver_type1,p1,xminmax1,c,dt1,custom1,F1,bound_type1);            
            tot_time_solver1 = tot_time_solver1 + toc(tstart);
            
            tstart = tic;
            %% n-1/2 |+| n -> n+1/2
            % SOLVE STATE 2
            [p2,n2] = solverWE1D_o(n2,solver_type2,p2,xminmax2,c,dt2,custom2,F2,bound_type2);
            tot_time_solver2 = tot_time_solver2 + toc(tstart);
            
            tstart = tic;
            % INTERFACE handling, compensating using (ONLY) PREVIOUS time step n
            % sim 1: n+1   (prev is n)
            % sim 2: n+1/2 (prev is n)
            % NOTE: compensate for interface pressures BEFORE temporal interpolation           
            [p1_next,p2_next] = solverWE1DInterface_o(p1,p2,...\
                x1d1,x1d2,dx1,dx2,dt1,dt2,c,x1d1_dx2,x1d2_dx1,order);
            p1(end-order/2+1:end,1) = p1_next;
            p2(1:order/2,1) = p2_next;
            tot_time_interface = tot_time_interface + toc(tstart);
            
            %% n |+| n+1/2 -> n+1
            tstart = tic;
            [p2,n2] = solverWE1D_o(n2,solver_type2,p2,xminmax2,c,dt2,custom2,F2,bound_type2);
            tot_time_solver2 = tot_time_solver2 + toc(tstart);
            
            tstart = tic;
            % create prev n+1/2 for sim 1 by interpolation            
            p1_dt2 = interp1(x,p1',xx,'cubic')'; % Note: p_current is not used for interface handling
            %p1_dt2_spline = spline(x,p1,xx); % Note: p_current is not used for interface handling
            %tot_time_tspline = tot_time_tspline + toc(tstart);
            
            %tstart = tic;
            % INTERFACE handling, compensating using (ONLY) PREVIOUS time step n+1/2
            % sim 1: n+1 (prev is n+1/2)
            % sim 2: n+1 (prev is n+1/2)
            %[~,p2_next] = solverWE1DInterface_o(p1_dt2',p2, ...\
            %    x1d1(1:5),x1d2,dx1,dx2,dt2,dt2,c,xminmax1,xminmax2,order);
            [~,p2_next] = solverWE1DInterface_o(p1_dt2,p2, ...\
                x1d1,x1d2,dx1,dx2,dt2,dt2,c,x1d1_dx2,x1d2_dx1,order);
            p2(1:order/2,1) = p2_next;
            tot_time_interface = tot_time_interface + toc(tstart);            
        else
            %% n-1 |+| n -> n+1
            % SOLVE (FULL) STATE 1
            [p1,n1] = solverWE1D_o(n1,solver_type1,p1,xminmax1,c,dt1,custom1,F1,bound_type1);
            
            %% n-1 |+| n -> n+1
            % SOLVE STATE 2
            [p2,n2] = solverWE1D_o(n2,solver_type2,p2,xminmax2,c,dt2,custom2,F2,bound_type2);
            
            % INTERFACE handling            
            [p1_next,p2_next] = solverWE1DInterface_o(p1,p2,x1d1,x1d2,dx1,dx2,dt1,dt2,c,x1d1_dx2,x1d2_dx1,order);
            p1(end-order/2+1:end,1)= p1_next;
            p2(1:order/2,1) = p2_next;
        end
        
        p1_all(n,:) = p1(:,1);
        p2_all(n,:) = p2(:,1);
        
        if mod(n-1,1000) == 0
            fprintf('Calculating n=%i ...\n',n-1)
        end
    end
    
    tot_time = toc(tstart_tot);
    
    fprintf('TOTAL: %0.3f \n', tot_time)
    fprintf('solverWE1D 1: %0.3f/%0.3f \n', tot_time_solver1, tot_time_solver1/tot_time)
    fprintf('solverWE1D 2: %0.3f/%0.3f \n', tot_time_solver2, tot_time_solver2/tot_time)
    fprintf('solverWE1DInterface: %0.3f/%0.3f \n', tot_time_interface, tot_time_interface/tot_time)
    %fprintf('adjustTemporal: %0.3f \n', tot_time_tspline)
end