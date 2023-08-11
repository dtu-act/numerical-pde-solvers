function [state] = solverWE1D_SEM(state)
    use_rk4 = false;
    
    F = state.src.F;
    n = state.n;

    c = state.domain.c;
    dt = state.domain.dt;

    p_current = state.p_current;
    p_prev = state.p_prev;
    
    Mx = state.custom.Mx;
    Sx = state.custom.Sx;
    conn = state.custom.conn;
    F_filter = state.custom.F_filter;
    
    switch state.domain.boundary
        case BoundaryType.Neumann
            % ok
        case BoundaryType.Dirichlet
            error("Dirichlet cond. not implemented")
        otherwise
            error('Boundary condition not supported');
    end
        
    if use_rk4
        dxP1_ = @(t, p) dxP1(p);
        dxP2_ = @(t, p) dxP2(p, Mx, Sx, c);
        [p_current, p_next] = rk4_coupled_weq(dxP2_, dxP1_, p_prev', p_current', n*dt, dt);
        p_next = p_next' + F(n*dt);
        p_current = p_current';
    else
        %rhs = -c^2 * Sx * p_current';
        %p_next = 2*p_current - p_prev + dt^2 * ( Mx \ rhs )' + F(n*dt);
        rhs = -c^2 * (Mx \ Sx) * p_current';
        p_next = 2*p_current - p_prev + dt^2 * rhs' + F(n*dt);
        nodes1 = conn(1,:); % conn(end,:)
        nodes2 = conn(2,:); % conn(end-1,:)
        p_next(nodes1) = F_filter*p_next(nodes1)'; % filter 1st element
        p_next(nodes2) = F_filter*p_next(nodes2)'; % filter 2nd element
    end
    
    state = state.update(p_next, p_current);
end

function [dP1] = dxP1(p2)
    dP1 = p2;
end

function [dP2] = dxP2(p1, Mx, Sx, c)
    rhs = - c^2 * Sx * p1;
    dP2 = Mx \ rhs;
end

function [p1, p2] = rk4_coupled_weq(dxP2, dxP1, p1, p2, t, dt)
    hHalf = dt/2.0;
    oneSixth = 1.0/6.0;

    K1p2 = dt * dxP2(t,          p1); % calculation of velocity is dependent on pressure
    K1p1 = dt * dxP1(t,          p2); % calculation of pressure is dependent on velocity

    K2p2 = dt * dxP2(t + hHalf,  p1 + K1p1./2.0);
    K2p1 = dt * dxP1(t + hHalf,  p2 + K1p2./2.0);

    K3p2 = dt * dxP2(t + hHalf,  p1 + K2p1./2.0);
    K3p1 = dt * dxP1(t + hHalf,  p2 + K2p2./2.0);

    K4p2 = dt * dxP2(t + dt,      p1 + K3p1);
    K4p1 = dt * dxP1(t + dt,      p2 + K3p2);

    p2 = p2 + oneSixth * (K1p2 + 2*(K2p2 + K3p2) + K4p2);
    p1 = p1 + oneSixth * (K1p1 + 2*(K2p1 + K3p1) + K4p1);
end
