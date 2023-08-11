function saveSimulationData(path_dir, Pfreqs, XY, xfreqs, t, freq_min_max, ...
    Porder, NeX, NeY, bbox, xy0_rel, boundary_type, aparams, tag)
    
    xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);
    path = paths.simulationDataPath(path_dir,tag,freq_min_max,Porder,NeX,NeY,aparams,xy0,bbox,boundary_type);
    
    if ~exist(path_dir, 'dir')
        mkdir(path_dir)
    end
    
    fprintf('Writing to %s\n', path);
    save(path, 'Pfreqs', 'XY', 'xfreqs', 't')
end