function [C] = constructConnectivity1D(etov, p)
%
% from element to global nodal index depending on polynomial order
% out: [1 2 3 4 5; 5 6 7 8 9] (for p=4, 2 elements)
%
    N = length(etov);
    Mp = p+1;
    C = zeros(N, Mp);
    gidx = 1;
    
    for n=1:N
        for i=1:Mp
            C(n,i) = gidx;
            gidx = gidx + 1;
        end                
        gidx = gidx - 1;
    end
end