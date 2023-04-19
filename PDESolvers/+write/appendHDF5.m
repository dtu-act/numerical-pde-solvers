function appendHDF5(filename,p,n)
    start = [1,1,n];
    count = [size(p,1),size(p,2),1];
    h5write(filename, '/pressures', single(p), start, count);
end