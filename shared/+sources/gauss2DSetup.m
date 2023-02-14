function [source] = gauss2DSetup(c,sigma,xy0)
    x0 = xy0(1);
    y0 = xy0(2);

    %sigma = c/(pi*f/2);
    F_gauss = @(x,y) -(exp(-(x-x0).^2/(sigma^2)) .* exp(-(y-y0).^2/(sigma^2)));
    source = models.SourceModel(models.SourceType.Function, F_gauss, xy0);
end