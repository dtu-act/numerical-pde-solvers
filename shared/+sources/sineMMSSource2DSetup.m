function [F, usol] = sineMMSSource2DSetup(k)
    usol = @(xy) sin(pi*xy(:,1)).*sin(pi*xy(:,2));
    F = @(x,y) -2*pi^2*sin(pi*x).*sin(pi*y) + k^2*sin(pi*x).*sin(pi*y);
end