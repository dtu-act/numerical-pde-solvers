function findBoundaries(p, etov, fd_inner, fd_outer)
    tol = h0/1e3;
    nodesInner = find(abs(fd_inner(p))<tol);
    nodesOuter = find(abs(fd_outer(p))<tol);
    
    BCType = int8(not(etov));
   
    % Choose a map to plot
    MAP = nodesInner;
    
    % Show all vertex nodes and circle out map nodes in red
    figure
    plot(p(:,1),p(:,2),'k.', p(MAP,1), p(MAP,2),'ro');
    
    In=1; Out=2;
    
    BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,nodesInner,fd_inner,In);
    BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,nodesOuter,fd_outer,Out);
        
	mapI = ConstructNodesMap(BCType,c,P,In); 
    mapO = ConstructNodesMap(BCType,c,P,Out);
    
    
end

function BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,fd,BCcode) 
% Purpose: Store boundary information in the BCType array
% EToV      : Element-To-Vertice table
% VX, VY    : (x,y)-coordinates of mesh vertices
% BCType    : Table with types of faces for BC's
% fd        : handle to distance function defining boundary
% BCcode    : Integer for specific boundary type
%
%    From notes "The Spectral/hp-Finite Element Method for Partial
%    Differential Equations" by by Allan P. Engsig-Karup
%
    VNUM = [1 2;2 3;3 1]; % face orientations
    pxc = 0.5*(VX(EToV)+VX(EToV(:,[2 3 1])));
    pyc = 0.5*(VY(EToV)+VY(EToV(:,[2 3 1])));
    dc = abs(fd([pxc(:) pyc(:)])); % distances to boundaries from face centers tol = 5e-2; % tolerance
    idx = find(dc<tol);
    BCType(idx) = BCcode;
end

function [map] = ConstructNodesMap(BCType,c,P,BCcode)
    Nfaces = size(BCType,2); 
    Npf = P+1; 
    [elms,fids] = find(BCType == BCcode); 
    map = zeros(length(elms)*Npf,1); 
    count = 0;
    
    for n = elms
        for f = fids
            lidx = [1:Nfaces Nfaces+(f-1)*(Npf-2)+(1:Npf-2)]; % local face indexes 
            map(count+(1:Npf)) = c(n,lidx); % vertex nodes
            count = count + Npf;
        end
    end
    
    map = unique(map);
end