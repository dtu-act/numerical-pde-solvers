function [state1_out, state2_out] = solverWE1DInterface(state1_in,state2_in,scheme_order)
    %%% immutable function
    
    % same parameters as input, but will have updates current pressures
    state1_out = copy(state1_in);
    state2_out = copy(state2_in);
    
    % used for intermediate modifications
    state1_tmp = copy(state1_in);
    state2_tmp = copy(state2_in);
    
    dx1 = getDx(state1_tmp);
    dx2 = getDx(state2_tmp);
    
    N_iface = scheme_order/2 + 1;
    
    if abs(dx1-dx2) > eps('double')
        % p prev from uniform grid with avg. resolution from right domain
        % used ONLY as 'ghost node' for right neighbour domain
        x1d_dx2 = fliplr(state1_tmp.domain.xmax:-dx2:(state1_tmp.domain.xmax - dx2*(N_iface-1))); % flip to ensure cont. at right interface
        p1p_i = spline(state1_tmp.domain.x1d,state1_tmp.p_prev,x1d_dx2);
        
        % p prev from uniform grid with avg. resolution from left domain
        % used ONLY as 'ghost node' for left neighbour domain    
        x1d_dx1 = state2_tmp.domain.xmin:dx1:(state2_tmp.domain.xmin + dx1*(N_iface-1));
        p2p_i = spline(state2_tmp.domain.x1d,state2_tmp.p_prev,x1d_dx1);
    else
        p1p_i = state1_tmp.p_prev;
        p2p_i = state2_tmp.p_prev;
    end
    
    % NOTE: p2p_i has same resolution as domain 1, p1p_i has same resolution as domain 2
    switch scheme_order        
        case 2
            solverWE1D2ndOrderInterface_Left(state1_tmp, p2p_i, dx1);
            solverWE1D2ndOrderInterface_Right(state2_tmp, p1p_i, dx2);
        case 6
            solverWE1D6thOrderInterface_Left(state1_tmp, p2p_i, dx1);
            solverWE1D6thOrderInterface_Right(state2_tmp, p1p_i, dx2);
        otherwise
            error('scheme order not implemented')
    end
    
    % update original states (including history)
    state1_out.adjust_current(state1_tmp.p_current);
    state2_out.adjust_current(state2_tmp.p_current);
end

function solverWE1D2ndOrderInterface_Left(state1, p2p_i, dx1)
    factor1 = state1.domain.c^2*state1.domain.dt^2/dx1^2;
    p1c = state1.p_current;
    p1p = state1.p_prev;
    
    % remove interface pressure terms for domain
    p1c(end) = p1c(end) - factor1*p1p(end-1);       % [1, -2, 1] -> [1, -2, 0]
    
    % add interface pressure terms from neighbour domain
    p1c(end) = p1c(end) + factor1*p2p_i(1+1);       % [1, -2, 0] -> [1, -2, 1]
    
    % adjust simulation states
    state1.adjust(p1c, p1p, false);
end

function solverWE1D2ndOrderInterface_Right(state2, p1p_i, dx2)
    factor2 = state2.domain.c^2*state2.domain.dt^2/dx2^2;
    p2c = state2.p_current;
    p2p = state2.p_prev;
    
    % remove interface pressure terms for domain
    p2c(1) = p2c(1) - factor2*p2p(1+1);       % [1, -2, 1] -> [0, -2, 1]
    
    % add interface pressure terms from neighbour domain
    p2c(1) = p2c(1) + factor2*p1p_i(end-1);   % [0, -2, 1] -> [1, -2, 1]
    
    % adjust simulation states
    state2.adjust(p2c, p2p, false);
end

function solverWE1D6thOrderInterface_Left(state, p2p_i, dx1)
    factor1 = state.domain.c^2*state.domain.dt^2/(180*dx1^2);
    
    p1c = state.p_current;
    p1p = state.p_prev;

    % remove interface pressure terms
    p1c(end)   = p1c(end)   - factor1*(270*p1p(end-1) - 27*p1p(end-2) + 2*p1p(end-3));
    p1c(end-1) = p1c(end-1) - factor1*(-27*p1p(end-1) + 2*p1p(end-2));
    p1c(end-2) = p1c(end-2) - factor1*(2*p1p(end-1));

    % add interface pressure terms from neighbour domain
    p1c(end)   = p1c(end)   + factor1*(270*p2p_i(1+1) - 27*p2p_i(1+2) + 2*p2p_i(1+3));
    p1c(end-1) = p1c(end-1) + factor1*(-27*p2p_i(1+1) + 2*p2p_i(1+2));
    p1c(end-2) = p1c(end-2) + factor1*(2*p2p_i(1+1));
    
    % adjust simulation states
    state.adjust(p1c, p1p, false);
end

function solverWE1D6thOrderInterface_Right(state, p1p_i, dx2)
    factor2 = state.domain.c^2*state.domain.dt^2/(180*dx2^2);
    
    p2c = state.p_current;
    p2p = state.p_prev;
    
    % remove interface pressure terms
    p2c(1)   = p2c(1)   - factor2*(270*p2p(1+1) - 27*p2p(1+2) + 2*p2p(1+3));
    p2c(1+1) = p2c(1+1) - factor2*(-27*p2p(1+1) + 2*p2p(1+2));
    p2c(1+2) = p2c(1+2) - factor2*(2*p2p(1+1));

    % add interface pressure terms from neighbour domain
    p2c(1)   = p2c(1)   + factor2*(270*p1p_i(end-1) - 27*p1p_i(end-2) + 2*p1p_i(end-3));
    p2c(1+1) = p2c(1+1) + factor2*(-27*p1p_i(end-1) + 2*p1p_i(end-2));
    p2c(1+2) = p2c(1+2) + factor2*(2*p1p_i(end-1));
    
    % adjust simulation states
    state.adjust(p2c, p2p, false);    
end