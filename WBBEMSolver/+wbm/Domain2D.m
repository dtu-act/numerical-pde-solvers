function [domain2D] = Domain2D(lx,ly,gcoord,nodes_bound,origin,boundary_type,source)
%Domain2D Initialize an Domain2D struct

    if nargin == 6
        source = models.SourceModelNone();
    end
    
    domain2D = wbm.DomainCoupling2D(lx, ly,...
        gcoord, nodes_bound, [], origin, boundary_type, source);
end