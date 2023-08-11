function [p,n] = solverWE1D_FOURIER_o(n,p,xminmax,c,dt,nmodes,F,bound_type)
    l = xminmax(2) - xminmax(1);
    
    pc = p(:,1)';
    pp = p(:,2)';
    
    k = @(n) pi*(n/l);
    w = @(n) c*k(n);
    ks = 0:(nmodes*2);
    assert(rem(nmodes*2,1) == 0)
    
    ws = w(ks);
    ws(1) = 1;

    M_current = dct.dctmodes(pc, bound_type);
    M_prev = dct.dctmodes(pp, bound_type); % speed optimization: we could keep the mode from prev (not modified)

    F_mode = dct.dctmodes(F(n*dt), bound_type); %F = highpass(F,20,source_max_freq*2); % highpass if pulse has a DC component

    % note: the forcing term is injecting little energy because of the
    % division with the angular frequencies
    F_mode_angular = 2*F_mode./ws.^2;
    M_next = 2*M_current.*cos(ws*dt) - M_prev + F_mode_angular.*(1 - cos(ws*dt));

    % without incorporating oscillation terms dt needs fine resolution to not get unstable
    %M_next = 2*M_current - M_prev + (-ws.^2.*M_current + F_mode)*dt^2; 

    pc = dct.idctmodes(M_next, bound_type);
    
    p(:,2:end) = p(:,1:end-1);
    p(:,1) = pc; % only the current pressure value is changed
                                                   
    n = n+1;
end