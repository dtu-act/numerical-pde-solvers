function [scatter, gather] = scatterGatherSEMnD(Nk,Np,gidx,conn)
%
%    Implemented by Allan P. Engsig-Karup, apek@dtu.dk.
%
%    Example:
%        _, gather = scatterGatherSEMnD(Nk,Np,gidx,conn)
%        x1d = gather*X[:]
%        y1d = gather*Y[:]
%
    
    % Create Global-To-Local (SCATTER) operator from connectivity tables in any spatial dimension.
    ii = zeros(1,Nk*Np);
    jj = zeros(1,Nk*Np);
    ss =  ones(1,Nk*Np);
    count = 0;

    for k = 1 : Nk
        for i = 1 : Np
            count = count + 1;
            ii(count) = count;
            jj(count) = conn(k,i);
        end
    end

    scatter = sparse(ii,jj,ss);
    
    % Create Local-To-Global (GATHER) operator
    % to a full 3D grid
    gather = sparse(gidx,Nk*Np);
    count = 0;
    for k = 1 : Nk
        for i = 1 : Np
            count = count + 1;
            gather(conn(k,i),:) = 0;
            gather(conn(k,i),count) = 1;
        end
    end
end