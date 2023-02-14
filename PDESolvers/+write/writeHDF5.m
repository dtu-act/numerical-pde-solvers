function writeHDF5(grid,t,p_sols,conn,x0_srcs,accumulators,...
    c,c_phys,rho,sigma,fmax,ppw,boundary_type,dx,domain_minmax,filename)

    if isfile(filename)
        delete(filename)
    end
    
    if ~isempty(x0_srcs)
        h5create(filename, '/x0_srcs', size(x0_srcs'), 'Datatype', 'single')
        h5write(filename, '/x0_srcs', single(x0_srcs)')
    end
    
    total_memory = (64*length(grid(:)) + 32*length(t(:)) + 32*length(p_sols(:)) + ...
        32*length(conn(:)) + 32*length(x0_srcs(:)) + 32*length(accumulators(:)))/8/1024/1024;
    fprintf('Total memory: %0.1f mb\n', total_memory)
    
    h5create(filename, '/p', size(p_sols), 'Datatype', 'single')
    h5write(filename, '/p', single(p_sols))
        
    %h5create(filename, '/v', size(v_sols), 'Datatype', 'single')
    %h5write(filename, '/v', v_sols)

    h5create(filename, '/grid', size(grid'))
    h5write(filename, '/grid', grid')

    h5create(filename, '/conn', size(conn'), 'Datatype', 'uint16')
    h5write(filename, '/conn', uint16(conn)')
    
    h5create(filename, '/t', length(t), 'Datatype', 'single')
    h5write(filename, '/t', single(t))
    
    if ~isempty(accumulators)
        h5create(filename, '/accumulators', size(accumulators), 'Datatype', 'single')
        h5write(filename, '/accumulators', single(accumulators))
    end

    dt = t(2)-t(1);

    fprintf('dx: %0.10f\n', dx)
    fprintf('dt: %0.10f\n', dt)
    
    h5writeatt(filename,'/','dt',dt)
    h5writeatt(filename,'/','dx',dx)
    h5writeatt(filename,'/','c',c)
    h5writeatt(filename,'/','c_phys',c_phys)
    h5writeatt(filename,'/','rho',rho)
    h5writeatt(filename,'/','sigma0',sigma)
    h5writeatt(filename,'/','fmax',fmax)
    h5writeatt(filename,'/','tmax',max(t,[],'all'))
    h5writeatt(filename,'/','ppw',ppw)
    h5writeatt(filename,'/','domain_minmax',domain_minmax)
    h5writeatt(filename,'/','boundary_type',boundary_type)
    
    %h5disp(filename);
end