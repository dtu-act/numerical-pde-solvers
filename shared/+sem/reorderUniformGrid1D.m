function [vxr, etovr] = reorderUniformGrid1D(vx, etov)
%
% vx:      assumes N uniformly distributed (but not necessarily ordered)
% etov:    assumes 1 to N-1 matching indexes in vx
% out:     Tuple{VX, EToV}
%          1: sorted element vertex node coordinates
%          2: sorted element to vertex table
%           

    vxr = zeros(1,length(vx));
    etovr = zeros(length(etov),2);
    
    for i=1:length(etov)
        [i0,~] = etov(i);
        vxr(i) = vx(i0);
        
        etovr(i,:) = [i,i+1];
    end
    
    [~,i1] = etov(end);
    vxr(end) = vx(i1);

    %return vxr, etovr
end