function [p1_next,p2_next] = solverWE1DInterface_o(p1,p2,x1d1,x1d2,dx1,dx2,dt1,dt2,c,x1d_dx2,x1d_dx1,scheme_order)
    %%% _o: optimized for speed
    
    p1c = p1(:,1);
    p1p = p1(:,2);
    p2c = p2(:,1);
    p2p = p2(:,2);
        
    if abs(dx1-dx2) > eps('double')
        % p prev from uniform grid with avg. resolution from right domain
        % used ONLY as 'ghost node' for right neighbour domain        
        p1p_i = spline(x1d1,p1p,x1d_dx2);
        %p1p_i = interp1(x1d1,p1p,x1d_dx2,'cubic');
        
        % p prev from uniform grid with avg. resolution from left domain
        % used ONLY as 'ghost node' for left neighbour domain            
        p2p_i = spline(x1d2,p2p,x1d_dx1);
    else
        p1p_i = p1p;
        p2p_i = p2p;
    end
    
    % NOTE: p2p_i has same resolution as domain 1, p1p_i has same resolution as domain 2
    switch scheme_order
        case 2
            [p1_next] = solverWE1D2ndOrderInterface_Left(p1c,p1p,p2p_i,dx1,dt1,c);
            [p2_next] = solverWE1D2ndOrderInterface_Right(p2c,p2p,p1p_i,dx2,dt2,c);
        case 6
            [p1_next] = solverWE1D6thOrderInterface_Left(p1c,p1p,p2p_i,dx1,dt1,c);
            [p2_next] = solverWE1D6thOrderInterface_Right(p2c,p2p,p1p_i,dx2,dt2,c);
        otherwise
            error('scheme order not implemented')
    end
    
    %p1(:,1) = p1c;
    %p2(:,1) = p2c;
end

function [p1c, p1p] = solverWE1D2ndOrderInterface_Left(p1c, p1p, p2p_i, dx1, dt1, c)
    factor1 = c^2*dt1^2/dx1^2;
    
    % remove interface pressure terms for domain
    p1c(end) = p1c(end) - factor1*p1p(end-1);       % [1, -2, 1] -> [1, -2, 0]
    
    % add interface pressure terms from neighbour domain
    p1c(end) = p1c(end) + factor1*p2p_i(1+1);       % [1, -2, 0] -> [1, -2, 1]
end

function [p2c, p2p] = solverWE1D2ndOrderInterface_Right(p2c, p2p, p1p_i, dx2, dt2, c)
    factor2 = c^2*dt2^2/dx2^2;
    
    % remove interface pressure terms for domain
    p2c(1) = p2c(1) - factor2*p2p(1+1);       % [1, -2, 1] -> [0, -2, 1]
    
    % add interface pressure terms from neighbour domain
    p2c(1) = p2c(1) + factor2*p1p_i(end-1);   % [0, -2, 1] -> [1, -2, 1]
end

function p1c = solverWE1D6thOrderInterface_Left(p1c, p1p, p2p_i, dx1, dt1, c)
    factor1 = c^2*dt1^2/(180*dx1^2);

    % remove interface pressure terms
    p1c_end   = p1c(end)   - factor1*(270*p1p(end-1) - 27*p1p(end-2) + 2*p1p(end-3));
    p1c_end_m1 = p1c(end-1) - factor1*(-27*p1p(end-1) + 2*p1p(end-2));
    p1c_end_m2 = p1c(end-2) - factor1*(2*p1p(end-1));

    % add interface pressure terms from neighbour domain
    p1c_end   = p1c_end   + factor1*(270*p2p_i(1+1) - 27*p2p_i(1+2) + 2*p2p_i(1+3));
    p1c_end_m1 = p1c_end_m1 + factor1*(-27*p2p_i(1+1) + 2*p2p_i(1+2));
    p1c_end_m2 = p1c_end_m2 + factor1*(2*p2p_i(1+1));
    
    p1c = [p1c_end_m2, p1c_end_m1, p1c_end];
end

function p2c = solverWE1D6thOrderInterface_Right(p2c, p2p, p1p_i, dx2, dt2, c)
    factor2 = c^2*dt2^2/(180*dx2^2);
    
    % remove interface pressure terms
    p2c_1 = p2c(1)   - factor2*(270*p2p(1+1) - 27*p2p(1+2) + 2*p2p(1+3));
    p2c_2 = p2c(1+1) - factor2*(-27*p2p(1+1) + 2*p2p(1+2));
    p2c_3 = p2c(1+2) - factor2*(2*p2p(1+1));

    % add interface pressure terms from neighbour domain
    p2c_1   = p2c_1   + factor2*(270*p1p_i(end-1) - 27*p1p_i(end-2) + 2*p1p_i(end-3));
    p2c_2 = p2c_2 + factor2*(-27*p1p_i(end-1) + 2*p1p_i(end-2));
    p2c_3 = p2c_3 + factor2*(2*p1p_i(end-1));
    
    p2c = [p2c_1,p2c_2,p2c_3];
end