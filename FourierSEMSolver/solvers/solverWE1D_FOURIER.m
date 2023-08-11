function [state] = solverWE1D_FOURIER(state)
    n = state.n;

    dt = state.domain.dt;
    c = state.domain.c;
    l = state.domain.xmax - state.domain.xmin;
    nmodes = state.custom.nmodes;

    boundcond = state.domain.boundary;

    k = @(n) pi*(n/l);
    w = @(n) c*k(n);
    ks = 0:(nmodes*2);
    assert(rem(nmodes*2,1) == 0)
    
    ws = w(ks);
    ws(1) = 1;

    M_current = dct.dctmodes(state.p_current, boundcond);
    M_prev = dct.dctmodes(state.p_prev, boundcond); % speed optimization: we could keep the mode from prev (not modified)

    F = state.src.F(n*dt);
    F_mode = dct.dctmodes(F, boundcond); %F = highpass(F,20,source_max_freq*2); % highpass if pulse has a DC component

    % note: the forcing term is injecting little energy because of the
    % division with the angular frequencies
    F_mode_angular = 2*F_mode./ws.^2;
    M_next = 2*M_current.*cos(ws*dt) - M_prev + F_mode_angular.*(1 - cos(ws*dt));

    % without incorporating oscillation terms dt needs fine resolution to not get unstable
    %M_next = 2*M_current - M_prev + (-ws.^2.*M_current + F_mode)*dt^2; 

    p_current = dct.idctmodes(M_next, boundcond);
    p_prev = state.p_current;

    state = state.update(p_current, p_prev);
end