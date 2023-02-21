function [p,vx,vy] = solveCoupledWaveEquation2D(z0,T,rho,N,M2D,Sx,Sy,ts,bc_update_fn,accs)
    if nargin == 9
        accs = [];
    end

    dzdt_ = @(t,z,accs) dzdt(t,z,M2D,Sx,Sy,rho,N,bc_update_fn,accs);
    z = sem.rk4_coupled_weq2(dzdt_, z0, T, ts, accs);

    vx = z(1:N,:);
    vy = z(N+1:2*N,:);
    p = z(2*N+1:3*N,:);
end

function [dz,rhs_accs] = dzdt(t,z,M2D,Sx,Sy,rho,N,bc_update_fn,accs)
    vx = z(1:N);       % x component of velocity field
    vy = z(N+1:2*N);   % y component of velocity field
    p = z(2*N+1:3*N);  % wave pressure    
        
    % Main system rhs
    rhs_vx = -1/rho*Sx*p;
    rhs_vy = -1/rho*Sy*p;
    
    rhs_accs = [];
    if isempty(accs)
        rhs_p = bc_update_fn(p,vx,vy);
    else
        [rhs_p, rhs_accs] = bc_update_fn(p,vx,vy,accs); % frequency dependent
    end
    
    dvx = M2D\rhs_vx;
    dvy = M2D\rhs_vy;
    dp  = M2D\rhs_p;    
    
    dz = [dvx;dvy;dp];
end