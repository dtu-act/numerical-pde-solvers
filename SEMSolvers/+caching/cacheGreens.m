function [P,XY] = cacheGreens(Porder,NeX,NeY,XY,...
    path_dir, aparams,xy0_rel,bbox,boundary_type,Nw_greens)

    xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);
    path = paths.refDataPath(path_dir,'GREENS',Porder,NeX,NeY,aparams,xy0,bbox,boundary_type,Nw_greens);
    P = caching.readWriteGreens(path,bbox,boundary_type,XY,xy0,aparams,Nw_greens);
end