function [XY,etov,bdetect_f] = mesh2D(bbox, NeX, NeY)
    xmin = bbox(1,1);
    ymin = bbox(1,2);
    xmax = bbox(2,1);
    ymax = bbox(2,2);

    [X, Y] = meshgrid(linspace(xmin,xmax,NeX + 1), linspace(ymin,ymax,NeY + 1));
    etov = meshing.createEToV(NeX,NeY);

    bdetect_f = meshing.createBoundRectDetectFunc(bbox, 1e-10);

    XY = [X(:), Y(:)];
end