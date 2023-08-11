function [dx, Nk] = calcResolutionEvenElements(l,fmax,ppw,P,c)
    dx = c/(ppw*fmax);
    dx_elem = dx*P;
    Nk = ceil(l/dx_elem);
    dx = (l/Nk)/P;
end