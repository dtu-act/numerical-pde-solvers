function source = gaussianSourceIC2D(sigma_x)
    source = @(x,y,xy0) exp(-((x-xy0(1)).^2/(sigma_x^2) + (y-xy0(2)).^2/(sigma_x^2)));
end