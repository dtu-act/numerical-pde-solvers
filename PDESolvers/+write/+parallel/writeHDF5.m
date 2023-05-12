function writeHDF5(filename,p,up,x0)
    fprintf('Writing to: %s\n',  filename)
%     if isfile(filename)
%         delete(filename)
%     end
    h5create(filename, '/pressures', size(p), 'Datatype', 'single')
    h5write(filename, '/pressures', single(p))

    h5create(filename, '/upressures', size(up), 'Datatype', 'single')
    h5write(filename, '/upressures', single(up))

    if ~isempty(x0)
        % empty for non-Gaussian IC types (e.g. GRF)
        h5create(filename, '/x0_srcs', size(x0'), 'Datatype', 'single')
        h5write(filename, '/x0_srcs', single(x0'))
    end
end