function initHDF5(filename,mesh,umesh,umesh_shape,p,up,tvec,conn,x0_srcs,...
    c,c_phys,rho,sigma,fmax,boundary_type,dx,dx_u,dx_src,domain_minmax)

    if isfile(filename)
        delete(filename)
    end
    
    p = permute(p,[2,3,1]);
    up = permute(up,[2,1]);

    h5create(filename, '/mesh', size(mesh'))
    h5write(filename, '/mesh', mesh')
    
    h5create(filename, '/umesh', size(umesh'))
    h5write(filename, '/umesh', umesh')
    h5writeatt(filename,'/umesh','umesh_shape',umesh_shape)
    
    total_memory = (64*length(mesh(:)) + 32*length(tvec(:)) + 32*length(p(:)) + ...
        32*length(conn(:)) + 32*length(x0_srcs(:)))/8/1024/1024;
    fprintf('Total memory: %0.1f mb\n', total_memory)
    
    h5create(filename, '/pressures', size(p), 'Datatype', 'single', 'Chunksize',[size(p,1) size(p,2) 1])
    h5create(filename, '/upressures', size(up), 'Datatype', 'single')
    h5write(filename, '/upressures', single(up))

    h5create(filename, '/conn', size(conn'), 'Datatype', 'uint16')
    h5write(filename, '/conn', uint16(conn)')
    
    h5create(filename, '/t', length(tvec), 'Datatype', 'single')
    h5write(filename, '/t', single(tvec))
    
    if ~isempty(x0_srcs)
        % empty for non-Gaussian IC types (e.g. GRF)
        h5create(filename, '/x0_srcs', size(x0_srcs'), 'Datatype', 'single')
        h5write(filename, '/x0_srcs', single(x0_srcs'))
    end

    dt = tvec(2)-tvec(1);

    fprintf('dx: %0.10f\n', dx)
    fprintf('dt: %0.10f\n', dt)
    
    h5writeatt(filename,'/upressures/','dx',dx_u)
    h5writeatt(filename,'/pressures/','dt',dt)
    h5writeatt(filename,'/pressures/','dx',dx)
    h5writeatt(filename,'/pressures/','dx_src',dx_src)
    h5writeatt(filename,'/pressures/','c',c)
    h5writeatt(filename,'/pressures/','c_phys',c_phys)
    h5writeatt(filename,'/pressures/','rho',rho)
    h5writeatt(filename,'/pressures/','sigma0',sigma)
    h5writeatt(filename,'/pressures/','fmax',fmax)
    h5writeatt(filename,'/pressures/','tmax',max(tvec,[],'all'))
    h5writeatt(filename,'/mesh/','domain_minmax',domain_minmax)
    h5writeatt(filename,'/mesh/','boundary_type',boundary_type)
    
    %h5disp(filename);
end