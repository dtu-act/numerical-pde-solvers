function [x1d, custom, dxMin, S, M, Dr] = setupSEMCustom(xminmax, dx, P_order, iface)        
    xmin = xminmax(1);
    xmax = xminmax(2);
    
    if isempty(iface)
        iface = SEMUniformInterface(0,0,1,false,"");        
    end
    
    [x1d,rs,rs_edge,nodes,etov,conn] = setupSEMMeshPrefine(xmin, xmax, P_order, dx, iface);    
    [Mx, Sx, V, dxMin, S, M, Dr] = assembly1D2ndOrderDerivativePrefine(x1d, rs, rs_edge, nodes, etov, conn, iface);
    
    %[x1d,rs,nodes,etov,conn] = setupSEMMesh(xmin, xmax, P_order, dx);
    %[Mx, Sx, V, dxMin, S, M, Dr] = assembly1D2ndOrderDerivative(x1d, P_order, rs, nodes, etov, conn);

    F = filter1D(V,P_order,0.9,4);

    custom = CustomSEM(Mx,Sx,V,F,conn,etov,P_order,iface);    
    x1d = x1d';
end

function [x1d,rs,nodes,etov,conn] = setupSEMMesh(xmin, xmax, P_order, dx)
    l = xmax - xmin;
    Nk = l / (dx*P_order);
    if Nk ~= round(Nk)
        disp('WARNING: Nk is rounded')
        assert( round(Nk) == Nk) % even number of elements
        Nk = round(Nk);
    else
        Nk = int32(Nk);
    end
    
    [nodes, etov, ~] = sem.setupUniformGrid1D(xmin, xmax, Nk+1);
    
    % local nodes as Legendre-Gauss-Lobatto (LGL) (p. 63)
    rs = spectral.JacobiGL(0, 0, P_order);
    %rs = linspace(-1,1,P_order+1)';
    
    % mesh1D requires sorted input: (nodes, etov) = reorderUniformGrid1D(nodes, etov)
    conn = sem.constructConnectivity1D(etov, P_order);

    % insert additional elements to the mesh, so that the LGL nodes gets inserted
    % into each local element
    x1d =  meshing.mesh1D(nodes, rs);
end

function [x1d,rs,rs_edge,vx,etov,conn] = setupSEMMeshPrefine(xmin, xmax, P_order, dx, iface)        
    l = xmax - xmin;
    
    dx_e = iface.dx;    
    P_order_e = iface.P_order;
    Nk_e = iface.Nk;
    
    % construct the vertice scaffold
    Nk = (l - dx_e*Nk_e*P_order_e) / (dx*P_order);
    

    if abs(Nk - round(Nk)) > eps('single')
        disp('WARNING: Nk is rounded')
        assert( round(Nk) == Nk) % even number of elements
        Nk = round(Nk);
    else
        Nk = int32(Nk);
    end
    
    if iface.loc == InterfaceLocation1D.LEFT
        vx_inner = linspace(xmin + dx_e*P_order_e*Nk_e, xmax, Nk+1);    
        vx_begin = linspace(xmin, xmin+dx_e*P_order_e*(Nk_e-1), Nk_e);        
        vx = [vx_begin,vx_inner];
    else
        vx_inner = linspace(xmin, xmax - dx_e*P_order_e*Nk_e, Nk+1);    
        vx_end = linspace(xmax - dx_e*P_order_e*(Nk_e-1), xmax, Nk_e);        
        vx = [vx_inner,vx_end];
    end    

    etov = sem.setupEToV1D(vx);
    
    % local nodes as Legendre-Gauss-Lobatto (LGL) (p. 63)
    rs = spectral.JacobiGL(0, 0, P_order);
    %rs = linspace(-1,1,P_order+1)';
    if iface.uniform_distr
        rs_edge = spectral.JacobiGL(0, 0, P_order_e);
    else
        rs_edge = linspace(-1,1,P_order_e+1)';
    end      
    
    % mesh1D requires sorted input
    conn = constructConnectivity1D(etov, P_order, P_order_e, iface);

    % insert additional elements to the mesh, so that the LGL nodes gets inserted
    % into each local element
    x1d =  mesh1D(vx, rs, rs_edge, iface);
end

function [vs1d] = mesh1D(vs, rs, rs_edge, iface)
%
% vs:      assumes N uniformly distributed (but not necessarily ordered)
% vs1d:    grid node coordinates including the nodes inside elements (p>1)
%
    N = length(vs);
    P = length(rs)-1;
    P_edge = length(rs_edge)-1;
    Nk_e = iface.Nk;
    
    vs1d = zeros((N-Nk_e-1)*P + Nk_e*P_edge,1) + 1;
    
    j = 1;
    for n=1:N-1
        offset = vs(n);
        h = vs(n+1) - vs(n);
        if iface.loc == InterfaceLocation1D.LEFT && n <= Nk_e ...\
                || iface.loc == InterfaceLocation1D.RIGHT && n >= N - Nk_e
            for i=1:P_edge+1                
                normR = (rs_edge(i)+1)/2.0; % r in (-1,1)
                vs1d(j + i - 1) = offset + normR*h;
            end
        else
            for i=1:P+1
                normR = (rs(i)+1)/2.0; % r in (-1,1)
                vs1d(j + i - 1) = offset + normR*h;                
            end            
        end
        j = j + i - 1;
    end
end

function [C] = constructConnectivity1D(etov,p1,p2,iface)
%
% from element to global nodal index depending on polynomial order
% out: [1 2 3 4 5; 5 6 7 8 9] (for p=4, 2 elements)
%
    N = length(etov);
    C = zeros(N, max(p1,p2)+1);
    gidx = 1;
    
    Nk_e = iface.Nk;
    
    for n=1:N
        if iface.loc == InterfaceLocation1D.LEFT && n <= Nk_e ...\
                || iface.loc == InterfaceLocation1D.RIGHT && n > N-Nk_e
            for i=1:p2+1
                C(n,i) = gidx;
                gidx = gidx + 1;
            end                
        else
            for i=1:p1+1
                C(n,i) = gidx;
                gidx = gidx + 1;
            end                
        end
        
        gidx = gidx - 1;
    end
end

function [Mx, Sx, V, dxMin, S, M, Dr] = assembly1D2ndOrderDerivative(x1d, P_order, rs, nodes, etov, conn)
    % Build reference element matrices
    V = spectral.Vandermonde1D(P_order, rs);
    Vr = spectral.GradVandermonde1D(P_order, rs);
    Dr = Vr / V;
    
    Np = length(rs);
    
    % general mass matrix
    M = inv(V*V');

    % general stiffness matrix
    S = Dr'*M*Dr;
    
    % step size between each node
    % NOTE: assumes that elements in EToV are sorted
    hs = nodes(etov(:,2)) - nodes(etov(:,1));
    
    [Mx, Sx] = sem.assembly.assembly1D(Np, hs, conn, M, S);
    
    % TIME STEPPING
    dxMin = min(x1d(2:end) - x1d(1:end-1)); % assume sorted
end

% see p. 130 in "Nodal Discontinuous Galerkin Methods", Jan Hesthaven 2013
function [F] = filter1D(V,N,cutoff,s)
    Nc = round(N*cutoff);
    invV = inv(V);

    filterdiag = ones(N+1,1);
    alpha = -log(eps);

    % initialize filter function
    for i=Nc:N
        filterdiag(i+1) = exp(-alpha*((i-Nc)/(N-Nc))^s);
    end

    F = V*diag(filterdiag)*invV;
end

function [Mx, Sx, V, dxMin, S, M, Dr] = assembly1D2ndOrderDerivativePrefine(x1d, rs, rs_edge, nodes, etov, conn, iface)
    P_order = length(rs)-1;
    P_order_e = length(rs_edge)-1;
    
    % Build reference element matrices
    V = spectral.Vandermonde1D(P_order, rs);
    Vr = spectral.GradVandermonde1D(P_order, rs);
    Dr = Vr / V;
    
    Ve = spectral.Vandermonde1D(P_order_e, rs_edge);
    Vre = spectral.GradVandermonde1D(P_order_e, rs_edge);
    Dre = Vre / Ve;
    
    Np = length(rs);
    
    % general mass matrix
    M = inv(V*V');
    Me = inv(Ve*Ve');
    
    % general stiffness matrix
    S = Dr'*M*Dr;
    Se = Dre'*Me*Dre;    
    
    % step size between each node
    % NOTE: assumes that elements in EToV are sorted
    hs = nodes(etov(:,2)) - nodes(etov(:,1));
    
    %[Mx, Sx] = sem.assembly.assembly1D(Nk, Np, hs, conn, M, S);
    [Mx, Sx] = assemblyPrefined1D(Np, hs, conn, M, S, Me, Se, iface);
    
    % TIME STEPPING
    dxMin = min(x1d(2:end) - x1d(1:end-1)); % assume sorted
end