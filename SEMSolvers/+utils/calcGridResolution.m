function [dx,NeX,NeY] = calcGridResolution(fmax, c, Porder, ppw, lx, ly)
    lambda = c/fmax;
    dx = lambda/ppw*Porder;
    NeX = max(ceil(lx/dx),2);
    NeY = max(ceil(ly/dx),2);
end