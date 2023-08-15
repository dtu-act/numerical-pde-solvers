function [p,n] = solverWE1D_SEM_o(n,p,c,dt,Mx,Sx,F,bound_type)
    import models.types.*
    
    use_rk4 = false;
    
    pc = p(:,1)';
    pp = p(:,2)';
    
    switch bound_type
        case BoundaryType.Neumann
            % ok
        case BoundaryType.Dirichlet
            error("Dirichlet cond. not implemented")
        otherwise
            error('Boundary condition not supported');
    end
        
    if use_rk4
        dxP1_ = @(t, p) dxP1(p);
        dxP2_ = @(t, p) dxP2(p, Mx, Sx, c, F(t)');
        [pp1, pc1] = rk4_coupled_weq(dxP2_, dxP1_, pp', pc', n, dt);
        pc = pc1';
        pp = pp1';
        
        p(:,2:end) = p(:,1:end-1);
        p(:,1) = pc;
        p(:,2) = pp;
    else
        %rhs = -c^2 * Sx * p_current';
        %p_next = 2*p_current - p_prev + dt^2 * ( Mx \ rhs )' + F(n);
        rhs = -c^2 * (Mx \ Sx) * pc';
        pc = 2*pc - pp + dt^2 * rhs' + F(n*dt);
        
        p(:,2:end) = p(:,1:end-1);
        p(:,1) = pc; % only the current pressure value is changed
    end
    
    n = n+1;
end

function [dP1] = dxP1(p2)
    dP1 = p2;
end

function [dP2] = dxP2(p1, Mx, Sx, c, Q)
    rhs = - c^2 * Sx * p1;
    dP2 = Mx \ rhs + Q;
end

function [p1, p2] = rk4_coupled_weq(dxP2, dxP1, p1, p2, n, dt)
    t = n*dt;
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
