function state = adjustResidual(state, scheme_order, loc)
    add_residual_part = false;        

    xx_in = state.domain.x1d; % use this when re-meshing to original mesh

    %  we assume SEM is uniformly distributed at interfaces
    if state.solver_type == SolverType.SEM
        if isempty(state.custom.interface) || state.custom.scheme_order > 2
            dx1 = state.domain.dx;
            xx_out = state.domain.xmin:dx1:state.domain.xmax;
            [state, ~] = state.remesh(xx_out, true); % history NOT remeshed                  
        else
            dx1 = state.custom.interface.dx;
        end
    else
        dx1 = state.domain.dx;
    end

    p1c = state.p_current;
    p1p = state.p_prev;
    
    dt1 = state.domain.dt;    
    c1 = state.domain.c;
        
    % remove interface pressure terms for each domain
    if state.domain.boundary == BoundaryType.Neumann
        switch scheme_order
            case 2
                factor = c1^2*dt1^2/dx1^2;
                
                if loc == InterfaceLocation1D.RIGHT
                    p1c(end) = p1c(end) - factor*p1p(end-1);   % [1, -2, 1] -> [1, -2, 0]
                else
                    p1c(1) = p1c(1) - factor*p1p(1+1);         % [1, -2, 1] -> [0, -2, 1]                
                end
                
                if add_residual_part
                    if loc == InterfaceLocation1D.RIGHT
                        p1c(end) = p1c(end) + factor*p1p(end);     % [1, -2, 0] -> [1, -1, 0]
                    else
                        p1c(1) = p1c(1) + factor*p1p(1);           % [0, -2, 1] -> [0, -1, 1]
                    end
                end
            case 6
                factor = c1^2*dt1^2/(180*dx1^2);
                
                if loc == InterfaceLocation1D.RIGHT
                    % remove residual part (ghost nodes)
                    p1c(end)   = p1c(end) - factor*(270*p1p(end-1) - 27*p1p(end-2) + 2*p1p(end-3));
                    p1c(end-1) = p1c(end-1) - factor*(-27*p1p(end-1) + 2*p1p(end-2));
                    p1c(end-2) = p1c(end-2) - factor*(2*p1p(end-1));
                else
                    p1c(1)   = p1c(1)   - factor*(270*p1p(1+1) - 27*p1p(1+2) + 2*p1p(1+3));
                    p1c(1+1) = p1c(1+1) - factor*(-27*p1p(1+1) + 2*p1p(1+2));
                    p1c(1+2) = p1c(1+2) - factor*(2*p1p(1+1));
                end
                
                if add_residual_part
                    if loc == InterfaceLocation1D.RIGHT
                        % add residual part inside domain
                        p1c(end)   = p1c(end)   + factor*(270*p1p(end) - 27*p1p(end-1) + 2*p1p(end-2));
                        p1c(end-1) = p1c(end-1) + factor*(-27*p1p(end) + 2*p1p(end-1));
                        p1c(end-2) = p1c(end-2) + factor*(2*p1p(end));
                    else
                        p1c(1)   = p1c(1)   + factor*(270*p1p(1) - 27*p1p(1+1) + 2*p1p(1+2));
                        p1c(1+1) = p1c(1+1) + factor*(-27*p1p(1) + 2*p1p(1+1));
                        p1c(1+2) = p1c(1+2) + factor*(2*p1p(1));
                    end
                end                
            otherwise
                error('scheme order not supported')
        end
    else
        error("Expected Neumann")
    end
    
    % add interface pressure terms for each domain    
    state = state.adjust(p1c, p1p);
    
    if state.solver_type == SolverType.SEM && isempty(state.custom.interface) || state.custom.scheme_order > 2
        state = state.remesh(xx_in, false);
    end  
end