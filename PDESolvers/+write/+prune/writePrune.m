function writePrune(todo)
    if boundary_type == "freq_indep"
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_xi%0.2f_%s%i.h5',...
            base_path,sim_id,fmax*c_phys,sigma0,c,xi,source_type,size(p_ics,1));
    elseif boundary_type == "freq_dep"
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_d%0.2f_%s%i.h5', ...
            base_path,sim_id,fmax*c_phys,sigma0,c,dmat,source_type,size(p_ics,1)); 
    else
        path_file = sprintf('%s/%s_2D_%0.2fHz_sigma%0.1f_c%i_%s%i.h5', ...
            base_path,sim_id,fmax*c_phys,sigma0,c,source_type,size(p_ics,1));
    end
    
    p_out = p_all; dx_out = dx;
    if ppw_x_out > -1
        [mesh,p_out,dx_out] = write.pruneSpatial(mesh,p_all,dx,ppw,ppw_x_out);        
    end
    
    dt_out = dt;
    if ppw_t_out > -1
        [tsteps_out,p_out,dt_out] = write.pruneTemporal(dt,fmax,tsteps_out,p_out,ppw_t_out);
    end
    
    % permute for .h5 format to be correct
    p_out_perm = permute(p_out,[2,3,1]);
    up_ics_perm = permute(up_ics,[2,3,1]);

    write.writeAllHDF5(mesh,umesh,umesh_shape,p_out_perm,up_ics_perm,tsteps_out,conn,x0_srcs,{},...
        c,c_phys,rho,sigma0,fmax,boundary_type,dx_out,[xminmax' yminmax'],path_file)
end