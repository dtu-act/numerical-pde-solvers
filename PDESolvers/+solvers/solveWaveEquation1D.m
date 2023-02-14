function [ps] = solveWaveEquation1D(T, c, pinit_1, pinit_2, Q, Mx, Sx, timesteps)
    dt = timesteps(2) - timesteps(1);
    use_rk4 = true;
    
    if use_rk4
        dxP1_ = @(t, p) dxP1(p, Q(t));
        dxP2_ = @(t, p) dxP2(p, Mx, Sx, c);    
        [ps, ~] = sem.rk4_coupled_weq(dxP2_, dxP1_, pinit_2, pinit_1, 0, T, timesteps);    
    else
        dxP_ = @(t, p) dxP(p, Mx, Sx, c);
        [ps] = sem.integrator_2nd_order(dxP_, pinit_1, pinit_2, dt, timesteps, Q);
    end    
end

function [dP1] = dxP1(p2, Q)
    dP1 = p2 + Q;
end

function [dP2] = dxP2(p1, Mx, Sx, c)
    rhs = - c^2 * Sx * p1;
    dP2 = Mx \ rhs;
end

function [dP] = dxP(p, Mx, Sx, c)
    rhs = - c^2 * Sx * p;
    dP = Mx \ rhs;
end