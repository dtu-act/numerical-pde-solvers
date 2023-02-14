function [Pref,XY] = cacheGreensK(aparams,path_dir,P,NeX,NeY,XY,xy0_rel,bbox,boundary_type,Nw_greens)
    path = paths.refDataPath(path_dir,'GREENS',P,NeX,NeY,aparams,xy0_rel,bbox,boundary_type,Nw_greens);        
    xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);
    Pref = caching.readWriteGreens(path,bbox,boundary_type,XY,xy0,aparams,Nw_greens);
end