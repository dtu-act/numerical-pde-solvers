function domain1d=Domain1D(x1d, dx, dt, c, rho, xminmax, boundary)
    domain1d = struct('dx',dx,'dt',dt,'c',c,'rho',rho,'x1d',x1d,...\
        'xmin',xminmax(1),'xmax',xminmax(2), 'l', xminmax(2)-xminmax(1), 'boundary',boundary);
end