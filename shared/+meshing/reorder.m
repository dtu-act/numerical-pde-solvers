function [etov] = reorder(etov,VX,VY)
% Purpose: Reorder elements to ensure counter-clockwise orientation 
%
%   From notes "The Spectral/hp-Finite Element Method for Partial
%   Differential Equations" by by Allan P. Engsig-Karup
%
    ax = VX(etov(:,1)); ay = VY(etov(:,1));
    bx = VX(etov(:,2)); by = VY(etov(:,2));
    cx = VX(etov(:,3)); cy = VY(etov(:,3));
    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    etov(i,:) = etov(i,[1 3 2]);
end