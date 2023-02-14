function [P,XY] = cacheSEM(Porder,NeX,NeY,path_dir,aparams,xy0_rel,bbox,boundary_type,only_element_nodes)
    xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);
    path = paths.refDataPath(path_dir,'SEM',Porder,NeX,NeY,aparams,xy0,bbox,boundary_type);
    [P,XY,conn] = caching.readWriteSEM(path,Porder,NeX,NeY,bbox,boundary_type,xy0,aparams);

    if only_element_nodes
        P = P(unique(conn(:,1:3)),:);
        XY = XY(unique(conn(:,1:3)),:);
    end
end