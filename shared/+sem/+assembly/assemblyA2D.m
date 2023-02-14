function [M, Sx, Sy] = assemblyA2D(Nk, Np, geofact, conn, Me, Dr, Ds)
%
% input:  Nk::Int, Np::Int, geofact::NamedTuple, conn::Matrix, M::Matrix, Dr::Matrix, Ds::Matrix
% output: Tuple{SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}}
%
% Nk   - number of elements
% Np   - number of nodes in each element = Nf + Nf*(P-1) (Nf=3 for triangles)
% Ntot - number of unique nodes - same as unique(etov)
%

    Ntot = length(unique(conn));
    M = zeros(Ntot, Ntot);
    Sx = zeros(Ntot, Ntot);
    Sy = zeros(Ntot, Ntot);
    
    rx = geofact.rx; sx = geofact.sx; ry = geofact.ry; sy = geofact.sy; J = geofact.J;
    
    tic;
    for k = 1:Nk
        % zero matrix with diagonal from Jacobian
        Jk = diag(J(:,k));
                
        % local mass matrix
        Mk = Jk * Me;
        
        % local stiffness matrix
        Dxk = diag(rx(:,k))*Dr + diag(sx(:,k))*Ds;
        Dyk = diag(ry(:,k))*Dr + diag(sy(:,k))*Ds;
        
        Skx = Jk * Dxk'*Me*Dxk;
        Sky = Jk * Dyk'*Me*Dyk;

%         for j = 1:Np
%             for i = 1:Np
%                 % get nodes making up k'th element
%                 indx_i = conn(k,i);
%                 indx_j = conn(k,j);
% 
%                 %if indx_j >= indx_i % exploit symmetry
%                 % update k'th element at index i,j
%                 %Mxy(:,indx_j) = sparse(indx_i,1,Mxy(indx_i,indx_j) + Mk(1:Np,j), Ntot, 1);
%                 M(indx_i,indx_j) = M(indx_i,indx_j) + Mk(i,j);
%                 Sx(indx_i,indx_j)  = Sx(indx_i,indx_j)  + Skx(i,j);
%                 Sy(indx_i,indx_j)  = Sy(indx_i,indx_j)  + Sky(i,j);
%                 %end
%             end            
%         end   
        
        for j = 1:Np
            % get nodes making up k'th element
            indx_i = conn(k,1:Np);
            indx_j = conn(k,j);

            M(indx_i,indx_j) = M(indx_i,indx_j) + Mk(1:Np,j);
            Sx(indx_i,indx_j)  = Sx(indx_i,indx_j)  + Skx(1:Np,j);
            Sy(indx_i,indx_j)  = Sy(indx_i,indx_j)  + Sky(1:Np,j);        
        end        
    end
    
    M = sparse(M);
    Sx = sparse(Sx);
    Sy = sparse(Sy);
    
    toc
end