function [EToE,EToF]= tiConnect2D(EToV,Nfaces)
% function [EToE,EToF]= tiConnect2D(EToV,Nfaces)
% Purpose: triangle face connect algorithm due to Toby Isaac
%          generalized by Allan P. Engsig-Karup

    K = size(EToV,1);
    Nnodes = max(max(EToV));

    % create list of all faces
    fnodes = [];
    for i = 1 : Nfaces
        if i==Nfaces
            fnodes = [fnodes; EToV(:,[i 1])];
        else
            fnodes = [fnodes; EToV(:,[i i+1])];
        end
    end

    fnodes = sort(fnodes,2)-1;

    % set up default element to element and Element to faces connectivity
    EToE = (1:K)'*ones(1,Nfaces);
    EToF = ones(K,1)*(1:Nfaces);

    % uniquely number each set of faces by their node numbers
    id = fnodes(:,1)*Nnodes + fnodes(:,2) + 1;
    spNodeToNode=[id, (1:Nfaces*K)', EToE(:), EToF(:)];

    % Now we sort by global face number.
    sorted=sortrows(spNodeToNode,1);

    % find matches in the sorted face list
    [indices,~]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

    % make links reflexive
    matchL = [sorted(indices,:)   ;sorted(indices+1,:)];
    matchR = [sorted(indices+1,:) ;sorted(indices,:)];

    % insert matches
    EToE(matchL(:,2)) = matchR(:,3);
    EToF(matchL(:,2)) = matchR(:,4);
end