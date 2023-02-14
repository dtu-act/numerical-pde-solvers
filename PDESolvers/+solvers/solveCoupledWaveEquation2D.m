function [p,vx,vy] = solveCoupledWaveEquation2D(z0,T,rho,N,M2D,Sx,Sy,ts,bc_update_fn)
    
    dzdt_ = @(t, z) dzdt(t,z,M2D,Sx,Sy,rho,N,bc_update_fn);
    z = sem.rk4_coupled_weq2(dzdt_, z0, T, ts);

    vx = z(1:N,:);
    vy = z(N+1:2*N,:);
    p = z(2*N+1:3*N,:);
end

function dz = dzdt(t,z,M2D,Sx,Sy,rho,N,bc_update_fn)
    vx = z(1:N);       % x component of velocity field
    vy = z(N+1:2*N);   % y component of velocity field
    p = z(2*N+1:3*N);  % wave pressure    
        
    % Main system rhs
    rhs_vx = -1/rho*Sx*p;
    rhs_vy = -1/rho*Sy*p;
    rhs_p = bc_update_fn(p,vx,vy);    
    
    dvx = M2D\rhs_vx;
    dvy = M2D\rhs_vy;
    dp  = M2D\rhs_p;    
    
    dz = [dvx;dvy;dp];
end