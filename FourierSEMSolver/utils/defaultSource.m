function source1d = defaultSource(fn_source, sourceFactor, domain, x0)
    nu = 1.3;
    %[source] = sources.gaussianSource1D(domain.c, fn_source);
    [source] = sources.gaussianSourceDiff1D(domain.c, fn_source*nu);
    F = @(t) source(domain.x1d,t,x0)*sourceFactor;    
    source1d = Source1D(F, x0);
end