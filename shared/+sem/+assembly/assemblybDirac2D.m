function [b] = assemblybDirac2D(Nk, Np, geofact, conn, M, x, y, source)
%
% input:  Nk::Integer, Np::Int, geofact::NamedTuple, conn::Matrix, M::Matrix, x::Vector, y::Vector, f::Function
% output: Vector
%
% Nk   - number of elements
% Np   - number of nodes in each element = Nf + Nf*(P-1) (Nf=3 for triangles)
% Ntot - number of unique nodes
%    
    Ntot = length(unique(conn));
    %Mxy = sparse(Ntot, Ntot);
    b = sparse(Ntot,1);
    
    x0 = source.r0(1);
    y0 = source.r0(2);
    Q = source.Q;
    
    elems_contributing = 0; % should be 6 of the point source is located on a grid point
    
    for k = 1:Nk        
        Jk = diag(geofact.J(:,k)); % zero matrix with diagonal from Jacobian
        %Mk = Jk * M; % local mass matrix
        Mk = 5/30*M; % local mass matrix
        did_contribute = false;
        
        for j = 1:Np
            jj = conn(k,j);
            xj = x(jj);
            yj = y(jj);
            
            if abs(xj - x0) < 1e-10 && abs(yj - y0) < 1e-10
                % we are inside elements contributing to the point source
                for i = 1:Np
                    ii = conn(k,i);
                    b(ii) = b(ii) + Mk(i,j)*Q;
                    
                    % DEBUGGING
                    did_contribute = true;                    
                    %bfull = full(b);
                    %fprintf('k = %i : ', k)
                    %fprintf('b(%i) = %i\n', ii, bfull(ii))
                    %Mxy(ii,jj) = Mxy(ii,jj) + Mk(i,j);
                end
            end         
        end
        
        if did_contribute
            elems_contributing = elems_contributing + 1;
        end
    end
    
    assert(elems_contributing == 6) % in case of point source        
    fprintf('elems contributing to point source: %i (==6 if structured grid and source is placed in an element node)\n', elems_contributing)    
end  