function [b] = assemblyb2D(Nk, Np, geofact, conn, M, x, y, F)
%
% Nk   - number of elements
% Np   - number of nodes in each element = Nf + Nf*(P-1) (Nf=3 for triangles)
% Ntot - number of unique nodes
%    
    Ntot = length(unique(conn));
    b = sparse(Ntot,1);
    
    for k = 1:Nk        
        Jk = diag(geofact.J(:,k)); % zero matrix with diagonal from Jacobian        
        Mk = Jk * M; % local mass matrix
        
        for j = 1:Np
            jj = conn(k,j);
            xj = x(jj);
            yj = y(jj);
            
            for i = 1:Np
                ii = conn(k,i);
                b(ii) = b(ii) + Mk(i,j)*F(xj,yj);
            end            
        end        
    end
end