function [Mx, Sx] = assembleWaveEquation1D(vx,etov,conn,rs)
    Porder = size(conn,2)-1;
    Np = length(rs);
    
    % Build reference element matrices
    V = spectral.Vandermonde1D(Porder, rs);
    Vr = spectral.GradVandermonde1D(Porder, rs);
    Dr = Vr / V;

    % general mass matrix
    M = inv(V*V');

    % general stiffness matrix
    S = Dr'*M*Dr;

    % step size between each node
    % NOTE: assumes that elements in EToV are sorted
    hs = vx(etov(:,2)) - vx(etov(:,1));
    
    [Mx, Sx] = sem.assembly.assembly1D(Np, hs, conn, M, S);