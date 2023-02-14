function path = refDataPath(path_dir,ref_type,P,NeX,NeY,aparam,xy0,bbox,boundary_type,Nw_greens)
    path_dir = paths.boundaryTypePath(path_dir,boundary_type);
    
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

    if nargin == 11
        path = sprintf('%s/%s%0.1fx%0.1f_P%i_Nex%i_Ney%i_k%0.2f_x%0.1fy%0.1f_alpha%0.1f_W%i.mat', ...
            path_dir,ref_type,lx,ly,P,NeX,NeY,aparam.k,xy0(1),xy0(2),alpha,Nw_greens);
    else
        path = sprintf('%s/%s%0.1fx%0.1f_P%i_Nex%i_Ney%i_k%0.2f_x%0.1fy%0.1f_alpha%0.1f.mat', ...
            path_dir,ref_type,lx,ly,P,NeX,NeY,aparam.k,xy0(1),xy0(2),alpha);
    end    
end