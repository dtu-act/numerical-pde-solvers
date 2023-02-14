function [L,b] = assemblyBounds1D(L,b,c0,c1,conn)
%
% output: Tuple{SparseMatrixCSC, Vector}
%

    N = size(conn,1);  % num elements
    Mp = size(conn,2); % num points per element
    
    % left boundary
    b(1) = c0;
    L(1,1:Mp) = [1, zeros(1,Mp-1)];
    
    % right boundary
    b(conn(N,Mp)) = c1;
    L(conn(N,Mp), conn(N,:)) = [zeros(1,Mp-1), 1];
end