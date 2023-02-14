function [vs1d] = mesh1D(vs, rs)
%
% vs:      assumes N uniformly distributed (but not necessarily ordered)
% etov:    assumes sorted from 1 to N-1 (matching indexes in vs)
% vs1d:    grid node coordinates including the nodes inside elements (p>1)
%
    N = length(vs);
    p = length(rs)-1;
            
    vs1d = zeros((N-1)*p+1,1);
    
    for n=1:N-1
        offset = vs(n);
        h = vs(n+1) - vs(n);
        for i=1:p+1
            normR = (rs(i)+1)/2.0; % r in (-1,1)
            vs1d((n-1)*p + i) = offset + normR*h;
        end
    end
end