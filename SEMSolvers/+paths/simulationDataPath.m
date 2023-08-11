function path = simulationDataPath(path_dir,tag,freq_min_max,P,NeX,NeY,aparam,xy0,bbox,boundary_type)
    path_dir = paths.boundaryTypePath(path_dir,boundary_type);
    
    fminmax_str = sprintf('fmin%0.1f_fmax%0.1f',freq_min_max(1),freq_min_max(2));
    
    xmin = bbox(1,1); xmax = bbox(2,1);
    ymin = bbox(1,2); ymax = bbox(2,2);
    
    lx =  xmax - xmin; ly =  ymax - ymin;
    
    if ~exist(path_dir, 'dir')
        mkdir(path_dir)
    end
    
    if boundary_type == models.BoundaryCondition.Impedance
        alpha = utils.ZtoAlpha(aparam.Z,aparam.c,aparam.rho);        
    else
        alpha = 1;
    end

    path = sprintf('%s/%s%0.1fx%0.1f_%s_P%i_Nex%i_Ney%i_x%0.1fy%0.1f_alpha%0.1f.mat', ...
        path_dir,tag,lx,ly,fminmax_str,P,NeX,NeY,xy0(1),xy0(2),alpha);  
end