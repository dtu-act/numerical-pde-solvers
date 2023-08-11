function run(iter, state, interface_order, do_adjust_residual)
    h = figure(1);    
    set(h,'Position',[10 10 1600 800])
    
    if do_adjust_residual
        L = state.domain.xmax - state.domain.min; 
        xx_out = state.domain.xmin:state.domain.dx:state.domain.xmax;
    end
    
    for n=1:iter        
        state = solverWE1D(state);
        if do_adjust_residual
            if false && state.solver_type == SolverType.SEM
                p_prev_init = state.p_prev;
                [state, cs1] = state.remesh(xx_out, true); % history NOT remeshed              
            end
            
            investigate_interp_before_iface = false;
            
            if investigate_interp_before_iface && n==230
                % investigate interpolation before interface
                cs = spline(state.domain.x1d, [0,state.p_current,0]);
                x = 0:0.005:0.2;
                p_interp = ppval(cs,x);
                figure(2); plot(state.domain.x1d, state.p_current, '-o'), hold on; plot(x, p_interp, '-')
            end
            
            state = adjustResidual(state, interface_order);
            
            if investigate_interp_before_iface && n==230
                % investigate interpolation before interface
                cs = spline(state.domain.x1d, [0,state.p_current,0]);
                x = 0:0.005:0.2;
                p_interp = ppval(cs,x);
                figure(3); plot(state.domain.x1d, state.p_current, '-o'), hold on; plot(x, p_interp, '-')
            end
        else
            p1 = state.p_current;
        end
        
        if false && state.solver_type == SolverType.SEM
            [state, ~] = state.remeshToInitial(true, p_prev_init); % history NOT remeshed
        end
        
        p1 = [0,state.p_current,0];
        x1d = [-state.domain.dx,state.domain.x1d,L + state.domain.dx];
            
        %plot_offset = 200+70;
        plot_offset = 1;
        nindx_zoom = 50;
        
        title_str = sprintf('%s, scheme order %i, interface order %i, n=%i', ...\
            state.solver_type, state.custom.scheme_order, interface_order, n);
        
         % include points +/- dx/2
        L1 = state.domain.x1d(end);
        
        x1d_1 = [-state.domain.dx, state.domain.x1d, L1 + state.domain.dx];        
        p1 = [0,state.p_current,0];
        
        if n >= 200 && mod(n,5) == 0
                plot(x1d_1,p1,'-ob')                
                
                if state.solver_type == SolverType.FOURIER
                    l1 = sprintf("%s, iface order %i",state.solver_type, interface_order);
                else
                    l1 = sprintf("%s P=%i, iface order %i",state.solver_type, state.custom.scheme_order, interface_order);
                end                
                
                legend(l1);
                title(title_str)
                xlabel('x')
                ylabel('p')
                ylim([-1,1])
                xlim([x1d_1(1),x1d_1(end)])
                drawnow
        end
            
%         if mod(n,5) == 0
%             h = figure(2);
%             set(h,'Position',[10 10 1600 800])
%             
%             x1d = [-state.domain.dx,state.domain.x1d,5+state.domain.dx];
%             subplot(2,1,1)
%             plot(x1d, [0,p_current,0], '-o');
%             ylim([-1,1])
%             xlim([x1d(1),x1d(end)])
%             
%             dxxP = dt^2 * ( Mx \ rhs )';            
%             subplot(2,1,2)
%             plot(state.domain.x1d, dxxP, '.');
%             
%             %pause(0.5)
%         end
%         
        %plotZoomed(n, x1d(plot_offset:end)', p1(plot_offset:end)', interface_order, nindx_zoom, ...\
        %    title_str, {sprintf("%s", state.solver_type)})
    end
end