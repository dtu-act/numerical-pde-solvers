function [ps,vs,laccs,raccs] = solveCoupledWaveEquation1D(T, c, rho, ps0, Q, A, Sx, timesteps, ...
    bc_update_fn, dtAccs)

    vs0 = ps0*0;

    % SOLVE
    dxV_ = @(t, p, laccs, raccs) dxV(p, A, Sx, rho);
    dxP_ = @(t, v, ps, laccs, raccs) dxP(v, A, Sx, rho, c, Q(t), bc_update_fn, ps, laccs, raccs);
    [ps, vs, laccs, raccs] = sem.rk4_coupled_weq(dxV_, dxP_, ps0, vs0, 0, T, timesteps, dtAccs);
end

function [dVs] = dxV(ps, A, Sx, rho)
    rhs = -1.0/rho * Sx*ps;
    dVs = A \ rhs;
end

function [dPs,laccs,raccs] = dxP(vs, A, Sx, rho, c, Q, bc_update_fn, ps, laccs, raccs)
    rhs = rho*c^2 * Sx'*vs;

    rhs(1) = bc_update_fn(rhs(1), ps(1), laccs);
    rhs(end) = bc_update_fn(rhs(end), ps(end), raccs);

    dPs = A \ rhs + Q;
end