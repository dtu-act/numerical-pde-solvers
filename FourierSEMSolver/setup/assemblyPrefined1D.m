function [Mx, Sx] = assemblyPrefined1D(Np, hs, connTable, M, S, Me, Se, iface)
    % Nk - number of elements
    % Np - number of nodes in each element = P+1
    
    Nk = size(connTable,1);
    Nk_e = iface.Nk;
    Np_edges = iface.P_order + 1;
    Nk_edges = iface.Nk;
    Nk_inner = Nk - Nk_edges;
    
    dim_edge = Nk_edges*Np_edges - Nk_edges;
    dim = Nk_inner*Np-Nk_inner+dim_edge+1;
    
    Mx = sparse(dim,dim);
    Sx = sparse(dim,dim);
    
    P = size(M,1)-1;
    Pe = size(Me,1)-1;
        
    for k = 1:Nk
        if iface.loc == InterfaceLocation1D.LEFT && k <= Nk_e ...\
                 || iface.loc == InterfaceLocation1D.RIGHT && k > Nk-Nk_e
           % local mass matrix
            Mk = hs(k)/2 * Me;
            % local stiffness matrix
            Sk = 2/hs(k) * Se;
            
            % get nodes making up k'th element
            idx = connTable(k,1:Pe+1);
        else        
            % local mass matrix
            Mk = hs(k)/2 * M;       
            % local stiffness matrix
            Sk = 2/hs(k) * S;
            
            % get nodes making up k'th element
            idx = connTable(k,1:P+1);
        end
        
        % update k'th element
        Mx(idx,idx) = Mx(idx,idx) + Mk; % TODO: slow implementation
        Sx(idx,idx) = Sx(idx,idx) + Sk;
    end
end