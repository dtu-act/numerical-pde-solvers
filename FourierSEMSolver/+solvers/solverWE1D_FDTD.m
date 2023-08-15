function [state] = solverWE1D_FDTD(state)
    switch state.custom.scheme_order
        case 2
            state = solverWE1D2ndOrder(state);
        case 6
            state = solverWE1D6thOrder(state);
        otherwise
            error('scheme order not implemented')
    end
end

function [state] = solverWE1D2ndOrder(state)
    import models.types.*
    
    n = state.n;

    c = state.domain.c;
    dt = state.domain.dt;
    dx = state.domain.dx;

    p_current = state.p_current;
    p_prev = state.p_prev;
    F = state.src.F(n*dt);

    switch state.domain.boundary
        case BoundaryType.Neumann
            p_current_m1 = [p_current(2),p_current(1:end-1)];   % [1, -2]
            p_current_p1 = [p_current(2:end),p_current(end-1)]; % [-2, 1]
        case BoundaryType.Dirichlet
            p_current_m1 = [0,0,p_current(2:end-1)];            % [0, -2]
            p_current_p1 = [p_current(2:end-1),0,0];            % [-2, 0]
        otherwise
            error('Boundary condition not supported');
    end
            
    factor = c^2*dt^2/dx^2;
    
    %p_next = (p_current_p1 + p_current_m1) - p_prev + F;
    p_next = 2*p_current - p_prev + factor*(p_current_p1 - 2*p_current + p_current_m1) + F;
    state = state.update(p_next, p_current);
end

function [state] = solverWE1D6thOrder(state)
    import models.types.*

    n = state.n;

    c = state.domain.c;
    dt = state.domain.dt;
    dx = state.domain.dx;

    p_current = state.p_current;
    p_prev = state.p_prev;
            
    switch state.domain.boundary
        case BoundaryType.Neumann
            p_current_m1 = [p_current(2),p_current(1:end-1)];
            p_current_m2 = [p_current(3),p_current(2),p_current(1:end-2)];
            p_current_m3 = [p_current(4),p_current(3),p_current(2),p_current(1:end-3)];
            p_current_p1 = [p_current(2:end),p_current(end-1)];
            p_current_p2 = [p_current(3:end),p_current(end-1),p_current(end-2)];
            p_current_p3 = [p_current(4:end),p_current(end-1),p_current(end-2),p_current(end-3)];
        case BoundaryType.Dirichlet
            p_current_m1 = [0,0,p_current(2:end-1)];
            p_current_m2 = [0,0,0,p_current(2:end-2)];
            p_current_m3 = [0,0,0,0,p_current(2:end-3)];
            p_current_p1 = [p_current(2:end-1),0,0];
            p_current_p2 = [p_current(3:end-1),0,0,0];
            p_current_p3 = [p_current(4:end-1),0,0,0,0];
        otherwise
            error('Boundary condition not supported');
    end
            
    factor = c^2*dt^2/(180*dx^2);
    
    p_next = factor*( .../
        2*p_current_m3 -27*p_current_m2 + 270*p_current_m1 .../
        - 490*p_current .../
        + 270*p_current_p1 - 27*p_current_p2 + 2*p_current_p3) + 2*p_current - p_prev + state.src.F(n*dt);
            
    state = state.update(p_next, p_current);
end
