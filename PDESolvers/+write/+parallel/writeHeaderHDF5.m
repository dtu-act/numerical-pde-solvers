function writeHeaderHDF5(filename,mesh,umesh,umesh_shape,tvec,conn,...
    c,c_phys,rho,sigma,fmax,boundary_type,dx,dx_u,dx_src,domain_minmax)

    if isfile(filename)
        delete(filename)
    end
    
    h5create(filename, '/mesh', size(mesh'))
    h5write(filename, '/mesh', mesh')
    
    h5create(filename, '/umesh', size(umesh'))
    h5write(filename, '/umesh', umesh')
    h5writeatt(filename,'/umesh','umesh_shape',umesh_shape)

    h5create(filename, '/conn', size(conn'), 'Datatype', 'uint16')
    h5write(filename, '/conn', uint16(conn)')
    
    h5create(filename, '/t', length(tvec), 'Datatype', 'single')
    h5write(filename, '/t', single(tvec))

    dt = tvec(2)-tvec(1);

    fprintf('dx: %0.10f\n', dx)
    fprintf('dt: %0.10f\n', dt)
    
    h5writeatt(filename,'/umesh/','dx',dx_u)
    h5writeatt(filename,'/mesh/','dt',dt)
    h5writeatt(filename,'/mesh/','dx',dx)
    h5writeatt(filename,'/mesh/','dx_src',dx_src)    
    h5writeatt(filename,'/mesh/','c',c)
    h5writeatt(filename,'/mesh/','c_phys',c_phys)
    h5writeatt(filename,'/mesh/','rho',rho)
    h5writeatt(filename,'/mesh/','sigma0',sigma)
    h5writeatt(filename,'/mesh/','fmax',fmax)
    h5writeatt(filename,'/mesh/','tmax',max(tvec,[],'all'))
    h5writeatt(filename,'/mesh/','domain_minmax',domain_minmax)
    h5writeatt(filename,'/mesh/','boundary_type',boundary_type)
    
    %h5disp(filename);
end