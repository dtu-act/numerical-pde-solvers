function [Mx, Sx, Hx] = assembly1D(Np, hs, connTable, M, S, H)
    Nk = size(connTable,1);

    % Nk - number of element
    % Np - number of nodes in each element = P+1    
    
    dim = Nk*Np-Nk+1;
    Mx = sparse(dim,dim);
    Sx = sparse(dim,dim);
    Hx = sparse(dim,dim);
    
    for k = 1:Nk
        % local mass matrix
        Mk = hs(k)/2 * M;       
        % local stiffness matrix
        Sk = 2/hs(k) * S;
        
        % get nodes making up k'th element
        idx = connTable(k,:);
        
        % update k'th element
        Mx(idx,idx) = Mx(idx,idx) + Mk; % TODO: slow implementation
        Sx(idx,idx) = Sx(idx,idx) + Sk; % TODO: slow implementation
        
        if nargin == 7
            %Hk = hs(k)/2 * H;
            Hk = H;
            Hx(idx,idx) = Hx(idx,idx) + Hk; % TODO: slow implementation
        end
    end
end