function [Pfreqs,XY,xfreqs,t] = loadSimulationData(path_dir, fminmax, ...
    Porder, NeX, NeY, bbox, xy0_rel, boundary_type, aparams, tag)

    xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);
    path = paths.simulationDataPath(path_dir,tag,fminmax,Porder,NeX,NeY,aparams,xy0,bbox,boundary_type);
    
    if ~exist(path_dir, 'dir')
        mkdir(path_dir)
    end
    
    fprintf('Loading from %s\n', path);
    data = load(path, 'Pfreqs', 'XY', 'xfreqs', 't');
    
    Pfreqs = data.Pfreqs;
    XY = data.XY;
    xfreqs = data.xfreqs;
    t = data.t;
end