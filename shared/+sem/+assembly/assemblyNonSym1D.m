function [Mx, Sx] = assemblyNonSym1D(Np, hs, connTable, M, S)
    Nk = size(connTable,1);
    
    % Nk - number of element
    % Np - number of nodes in each element = P+1
    
    dim = Nk*Np-Nk+1;
    Mx = sparse(dim,dim);
    Sx = sparse(dim,dim);
    
    for k = 1:Nk
        % local mass matrix
        Mk = hs(k)/2 * M;
        % local stiffness matrix
        Sk = S;
        
        % get nodes making up k'th element
        idx = connTable(k,:);
        
        % update k'th element
        Mx(idx,idx) = Mx(idx,idx) + Mk;
        Sx(idx,idx) = Sx(idx,idx) + Sk;
    end
end
