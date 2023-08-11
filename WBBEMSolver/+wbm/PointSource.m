function [pointSource2D] = PointSource(r0,Q)
%POINTSOURCE Initialize a PointSource struct

    dpqhat_fun = wbm.dpqhatPoint2D(r0(1),r0(2),Q); % derivative
    pqhat_fun = wbm.pqhatPoint2D(r0(1),r0(2),Q);   % pressure
    
    pointSource2D = wbm.SourceModel(models.SourceType.PointSource, ...
        r0, pqhat_fun, dpqhat_fun);
end