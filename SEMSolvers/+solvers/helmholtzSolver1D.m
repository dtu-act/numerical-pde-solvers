function [ps, x1d] = solveHelmholtz1D(L, k, P_order, c, fmax, ppw, source, bounds)
    c0 = bounds(1);
    c1 = bounds(2);

    xl = 0.0;
    xu = L;
   
    % SETUP THE MESH
    % NOTE: we add +1 to the number of nodes (dof), otherwise the implementation logic in assembly fails
    [nodes, etov, Ms] = sem.setupUniformGrid1D(c, fmax, P_order, ppw, xl, xu);

    % local nodes as Legendre-Gauss-Lobatto (LGL) (p. 63)
    rs = spectral.JacobiGL(0, 0, P_order);
    Np = length(rs);

    conn = meshing.constructConnectivity1D(etov, P_order);

    % insert additional elements to the mesh, so that the LGL nodes gets inserted
    % into each local element
    x1d =  meshing.mesh1D(nodes, rs);

    % ASSEMBLY
    % Build reference element matrices
    V = spectral.Vandermonde1D(P_order, rs);
    Vr = spectral.GradVandermonde1D(P_order, rs);
    Dr = Vr / V;

    % general mass matrix
    M = inv(V*V');

    % general stiffness matrix
    K = Dr'*M*Dr;

    H = spectral.Lagrange1D(rs,rs);
    
    % step size between each node
    % NOTE: assumes that elements in EToV are sorted
    hs = nodes(etov(:,2)) - nodes(etov(:,1));
    
    [Mk, Kk, Hk] = sem.assembly.assembly1D(Np, hs, conn, M, K, H);
    
    L = -Kk + k^2*Mk;
    
    switch source.type
        case models.SourceType.Function
            F  = source.F;
            b = Mk*F(x1d);
        case models.SourceType.PointSource
            F = dirac1D(source.r0,source.Q,1e-10);
            b = Hk*F(x1d);
        otherwise
            error('Source type not supported')
    end

    [Ltilde, btilde] = sem.assembly.assemblyBounds1D(L,b,c0,c1,conn);
    
    ps = Ltilde \ btilde;
end

function dirac_f = dirac1D(x0,Q,precision)
    function [vals] = dirac(x)
        vals = zeros(length(x),1);  
        indxs = find(abs(x-x0) < precision);
        vals(indxs) = Q;
    end

    dirac_f = @dirac;
end