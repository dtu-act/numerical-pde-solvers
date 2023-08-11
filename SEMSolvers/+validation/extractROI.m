function [Proi] = extractROI(P,XY,dist_f)
    d = dist_f(XY);
    indxs = find(d < 0);
    
    Proi = P;
    Proi(indxs) = 0;
end