function [Mxy] = assemblyMsub2D(Nk, Np, geofact, conn, M, x, y, x0, y0)
%
% input:  Nk::Int, Np::Int, geofact::NamedTuple, conn::Matrix, M::Matrix, Dr::Matrix, Ds::Matrix
% output: Tuple{SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}}
%
% Nk   - number of elements
% Np   - number of nodes in each element = Nf + Nf*(P-1) (Nf=3 for triangles)
% Ntot - number of unique nodes - same as unique(etov)
%

    Ntot = length(unique(conn));
    Mxy = sparse(Ntot, Ntot);
    
    J = geofact.J;

%     for k = 1:Nk        
%         Jk = diag(geofact.J(:,k)); % zero matrix with diagonal from Jacobian
%         Mk = Jk * M; % local mass matrix
% 
%         for j = 1:Np
%             jj = conn(k,j);
%             xj = x(jj);
%             yj = y(jj);
%             for i = 1:Np
%                 ii = conn(k,i);
%                 %b(ii) = b(ii) + Mk(i,j)*F(xj,yj);
%                 
%                 b(ii) = b(ii) + Mk(i,j)*F(xj,yj);
%             end            
%         end        
%     end
    elems_contributing = 0;
    
    for k = 1:Nk
        % zero matrix with diagonal from Jacobian
        Jk = diag(J(:,k));
                
        % local mass matrix
        Mk = Jk * M;
        
        delta_x = abs( x(conn(k,:)) - x0 );
        delta_y = abs( y(conn(k,:)) - y0 );
        
        if isempty( find(delta_x < 1e-10, 1 ) ) || isempty( find(delta_y < 1e-10, 1 ))
            % only collect contributions for elements where (x0,y0) is
            % included
            continue
        end
        
        elems_contributing = elems_contributing + 1;
        for j = 1:Np
            for i = 1:Np
                % get nodes making up k'th element
                ii = conn(k,i);
                jj = conn(k,j);
                
                fprintf('(x,y) = (%0.2f, %0.2f)\n', x(ii), y(jj))
                
                %if indx_j >= indx_i % exploit symmetry
                % update k'th element at index i,j
                Mxy(ii,jj) = Mxy(ii,jj) + Mk(i,j); % TODO: slow implementation
                %end
            end
        end
        fprintf('---- \n')
    end
    
    fprintf('elems_contributing: %i\n', elems_contributing)
end