function [dx, nmodes] = calcSpatial(l,fmax,ppw,c)
    lambda_min = c/fmax;
    
    refine_factor = (ppw/2);
    nmodes = ceil(l/lambda_min)*refine_factor;
    dx = l/(2*nmodes);
end