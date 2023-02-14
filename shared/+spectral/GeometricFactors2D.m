function geo_factors = GeometricFactors2D(x,y,Dr,Ds)

% Purpose  : Compute the metric elements for the local mappings of the elements
% Output   : struct with fields {rx | sx | ry | sy | J | xr | xs | yr | ys} 

    % Calculate geometric factors
    xr = Dr*x; 
    xs = Ds*x; 
    yr = Dr*y; 
    ys = Ds*y; 

    J = -xs.*yr + xr.*ys;

    rx = ys./J; 
    sx =-yr./J; 
    ry =-xs./J; 
    sy = xr./J;
          
    geo_factors = struct(...
        'xr', xr, 'xs', xs, 'yr', yr, 'ys', ys, ...
        'J', J, ...
        'rx', rx, 'sx', sx, 'ry', ry, 'sy', sy);
end


