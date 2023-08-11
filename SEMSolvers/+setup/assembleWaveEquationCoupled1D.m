function [A, Sx] = assembleWaveEquationCoupled1D(vx,etov,conn,rs)   
    Porder = size(conn,2)-1;
    
    % ASSEMBLY    
    % Build reference element matrices
    V = spectral.Vandermonde1D(Porder, rs);
    Vr = spectral.GradVandermonde1D(Porder, rs);
    Dr = Vr / V;
    
    % general mass matrix
    M = inv(V*V');
    
    % general stiffness matrix
    S = M*Dr;
    
    % step size between each node
    % NOTE: assumes that elements in EToV are sorted
    hs = vx(etov(:,2)) - vx(etov(:,1));    
    [A, Sx] = sem.assembly.assemblyNonSym1D(length(rs), hs, conn, M, S);
end