function dt = calculateDt(r,s,x,y,CFL,c,P,NODETOL)
    rLGL = spectral.JacobiGL(0,0,P);
    rmin = abs(rLGL(1)-rLGL(2));
    dtscale = min(setup.dtscale2D(r,s,x,y,NODETOL));    
    dt = min(dtscale)*rmin*CFL/c;

    assert(dt > 0)