function [Pdomains] = greensSpatialDomains2D(X, Y, domains, k, Nw_greens)
    if nargin == 4
        Nw_greens = 200;
    end
    
    Pdomains = cell(length(domains));
    for i=1:length(domains)
        domain = domains{i};
        Xi = X{i};
        Yi = Y{i};
        XYi = [Xi(:), Yi(:)];
        r0 = domain.source.r0;
        
        if domain.boundary_type == models.BoundaryCondition.Pressure
            [P, ~] = greens.greensDirichletUnstructuredSpatial2D(Nw_greens,Nw_greens,XYi,r0,k);
        else
            [P, ~] = greens.greensNeumannUnstructuredSpatial2D(Nw_greens,Nw_greens,XYi,r0,k);
        end
        
        Pdomains{i} = P;
    end
end