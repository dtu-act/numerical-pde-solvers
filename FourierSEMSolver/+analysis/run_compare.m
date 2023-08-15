function run_compare(iter, state1, state2, adjust, interface_order1, interface_order2, interface_loc)
    import solvers.*
    import models.types.*
    import plotting.*

    plot_zoomed = false;
    write_gif = true;
    plot_offset = 0; %900
    
    h = figure(1);
    set(h,'Position',[10 10 1600 800])    

    for n=1:iter
        switch adjust
            case "spatial"
                state1 = solverWE1D(state1);
                state2 = solverWE1D(state2);
                state2 = solverWE1D(state2);
            case "temporal"
                state1 = solverWE1D(state1);
                state1 = adjustResidual(state1, interface_order1, interface_loc);
                state2 = solverWE1D(state2);
                state2 = adjustResidual(state2, interface_order2, interface_loc);
                state2 = solverWE1D(state2);
                state2 = adjustResidual(state2, interface_order2, interface_loc);
            case "none"
                test_filter = true;

                state1 = solverWE1D(state1);
                if test_filter
                    x1d_uniform = state2.domain.x1d;
                    x1d = state1.domain.x1d;
                    dx1 = state1.domain.dx;
                    dt1 = state1.domain.dt;
                    c1 = state1.domain.c;
                    p1_current = adjustResidual_o(state1.p_current,state1.p_prev,x1d,dx1,dt1,c1,x1d_uniform,interface_order1,interface_loc);
                    state1 = state1.adjust_current(p1_current);
                else
                    state1 = adjustResidual(state1, interface_order1, interface_loc);
                end
                
                state2 = solverWE1D(state2);                       
                state2 = adjustResidual(state2, interface_order2, interface_loc);
            otherwise
                error("Compare type not implemented")
        end
        
        if state1.solver_type == SolverType.SEM && ~isempty(state1.custom.interface)
            dx1 = state1.custom.interface.dx;
        else
            dx1 = state1.domain.dx;
        end
        
        if state2.solver_type == SolverType.SEM && ~isempty(state2.custom.interface)
            dx2 = state2.custom.interface.dx;
        else
            dx2 = state2.domain.dx;
        end
        
        % include points +/- dx/2
        x1d_1 = [-dx1 + state1.domain.xmin, state1.domain.x1d];
        x1d_2 = [-dx2 + state2.domain.xmin, state2.domain.x1d];
        p1 = [0,state1.p_current];
        p2 = [0,state2.p_current];
        
        title_str = sprintf('n=%i', n);        
        
        if plot_zoomed
            legends = {sprintf("%s, iface order %i", state1.solver_type, interface_order1),sprintf("%s, iface order %i", state2.solver_type, interface_order2)};
            if n > plot_offset
                plotZoomed(n, [x1d_1', x1d_2'],  ...\
                              [p1(plot_offset:end)', p2(plot_offset:end)'], [interface_order1,interface_order2], ...\
                              nindx_zoom, plot_offset, title_str, legends)
               if write_gif
                   filename = '~/data/fm-sem/solvers.gif';
                   writeGif(h, filename, n==plot_offset);
               end
            end
        else
            if n >= plot_offset && mod(n-1,5) == 0
                subplot(2,1,1)
                plot(x1d_1,p1,'-o')
                hold on
                plot(x1d_2,p2,'-*')
                hold off
                
                if state1.solver_type == SolverType.FOURIER
                    l1 = sprintf("%s, iface order %i",state1.solver_type, interface_order1);
                else
                    l1 = sprintf("%s P=%i, iface order %i",state1.solver_type, state1.custom.scheme_order, interface_order1);
                end
                
                if state2.solver_type == SolverType.FOURIER
                    l2 = sprintf("%s, iface order %i", state2.solver_type, interface_order2);
                else
                    l2 = sprintf("%s P=%i, iface order %i",state2.solver_type, state2.custom.scheme_order, interface_order2);
                end     
                
                legend(l1,l2);
                title(title_str)
                xlabel('x')
                ylabel('p')
                ylim([-1,1])
                xlim([min(x1d_1(1),x1d_2(1)),max(x1d_1(end),x1d_2(end))])
                drawnow
            
                if state1.solver_type == SolverType.SEM
                    dxxP1 = (state1.custom.Mx \ state1.custom.Sx) * state1.p_current';
                else
                    p1_current = state1.p_current;
                    p1_current_m1 = [state1.p_current(2),state1.p_current(1:end-1)];   % [1, -2]
                    p1_current_p1 = [state1.p_current(2:end),state1.p_current(end-1)]; % [-2, 1]
                    dxxP1 = -1/state1.domain.dx^2*(p1_current_p1 - 2*p1_current + p1_current_m1)';
                end
                
                if state2.solver_type == SolverType.SEM
                    dxxP2 = (state2.custom.Mx \ state2.custom.Sx) * state2.p_current';
                else
                    p2_current = state2.p_current;
                    p2_current_m1 = [state2.p_current(2),state2.p_current(1:end-1)];   % [1, -2]
                    p2_current_p1 = [state2.p_current(2:end),state2.p_current(end-1)]; % [-2, 1]
                    dxxP2 = -1/state2.domain.dx^2*(p2_current_p1 - 2*p2_current + p2_current_m1)';
                end
                
                %dxxP1 = (state1.custom.Mx \ state1.custom.Sx) * state1.p_current';
                %dxxP2 = (state2.custom.Mx \ state2.custom.Sx) * state2.p_current';
                
                %dxxP1 = (state1.custom.Sx) * state1.p_current';
                %dxxP2 = (state2.custom.Sx) * state2.p_current';
                
                %dxxP1 = state1.domain.dt^2 * ( state1.custom.Mx \ rhs1 )';    
                %dxxP2 = state2.domain.dt^2 * ( state2.custom.Mx \ rhs2 )';    
                
                subplot(2,1,2)
                plot(state1.domain.x1d,dxxP1,'-o')
                hold on
                plot(state2.domain.x1d,dxxP2,'-o')
                hold off                
                legend(l1,l2);
                title('Laplacian term comparison for FDTD and SEM', 'fontsize', 16)
                
%                 if n == 921
%                     filename = sprintf('~/data/fm-sem/residual_removed_order1_%i_order2_%i_n%i.mat', ...\
%                         state1.custom.scheme_order,state2.custom.scheme_order,n);
%                     s1_type = state1.solver_type;
%                     s2_type = state2.solver_type;
%                     
%                     save(filename, 'x1d_1', 'x1d_2', 'p1', 'p2', 'dxxP1', 'dxxP2', 's1_type', 's2_type')
%                 end
                
                if write_gif
                    filename = '~/data/fm-sem/residual_removed.gif';
                    writeGif(h, filename, n==1);
                end
            end
        end        
    end
    filename = sprintf('~/data/fm-sem/residual_removed_order1_%i_order2_%i_n%i.mat', ...\
            state1.custom.scheme_order,state2.custom.scheme_order,n);
    s1_type = state1.solver_type;
    s2_type = state2.solver_type;
    
    save(filename, 'x1d_1', 'x1d_2', 'p1', 'p2', 'dxxP1', 'dxxP2', 's1_type', 's2_type')
end