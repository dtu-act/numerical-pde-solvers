function p1c = adjustResidual_o(p1c,p1p,x1d1,dx1,dt1,c,x1d1_uniform,scheme_order,loc)
    %%% _o: optimized for speed
    
    assert(abs(dx1 - (x1d1_uniform(2) - x1d1_uniform(1))) < eps('double'))

    p1c_i = spline(x1d1,p1c,x1d1_uniform);
    p1p_i = spline(x1d1,p1p,x1d1_uniform);
    
    switch scheme_order
        case 2
            if loc == InterfaceLocation1D.LEFT
                p1_uniform = solverWE1D2ndOrderToDirichlet_Left(p1c_i,p1p_i,dx1,dt1,c);
            else
                p1_uniform = solverWE1D2ndOrderToDirichlet_Right(p1c_i,p1p_i,dx1,dt1,c);               
            end
        case 6
            if loc == InterfaceLocation1D.LEFT
                p1_uniform = solverWE1D6thOrderToDirichlet_Left(p1c_i,p1p_i,dx1,dt1,c);
            else
                p1_uniform = solverWE1D6thOrderToDirichlet_Right(p1c_i,p1p_i,dx1,dt1,c);               
            end
        otherwise
            error('scheme order not implemented')
    end

    p1c(1) = p1_uniform(1);
end

function p1c = solverWE1D2ndOrderToDirichlet_Right(p1c, p1p, dx1, dt1, c)
    factor1 = c^2*dt1^2/dx1^2;
    
    % remove interface pressure terms for domain
    p1c(end) = p1c(end) - factor1*p1p(end-1);       % [1, -2, 1] -> [1, -2, 0]
end

function p2c = solverWE1D2ndOrderToDirichlet_Left(p2c, p2p, dx2, dt2, c)
    factor2 = c^2*dt2^2/dx2^2;
    
    % remove interface pressure terms for domain
    p2c(1) = p2c(1) - factor2*p2p(1+1);       % [1, -2, 1] -> [0, -2, 1]
end

function p1c = solverWE1D6thOrderToDirichlet_Right(p1c, p1p, dx1, dt1, c)
    factor1 = c^2*dt1^2/(180*dx1^2);

    % remove interface pressure terms
    p1c_end   = p1c(end)   - factor1*(270*p1p(end-1) - 27*p1p(end-2) + 2*p1p(end-3));
    p1c_end_m1 = p1c(end-1) - factor1*(-27*p1p(end-1) + 2*p1p(end-2));
    p1c_end_m2 = p1c(end-2) - factor1*(2*p1p(end-1));
    
    p1c = [p1c_end_m2, p1c_end_m1, p1c_end];
end

function p2c = solverWE1D6thOrderToDirichlet_Left(p2c, p2p, dx2, dt2, c)
    factor2 = c^2*dt2^2/(180*dx2^2);
    
    % remove interface pressure terms
    p2c_1 = p2c(1)   - factor2*(270*p2p(1+1) - 27*p2p(1+2) + 2*p2p(1+3));
    p2c_2 = p2c(1+1) - factor2*(-27*p2p(1+1) + 2*p2p(1+2));
    p2c_3 = p2c(1+2) - factor2*(2*p2p(1+1));
    
    p2c = [p2c_1,p2c_2,p2c_3];
end