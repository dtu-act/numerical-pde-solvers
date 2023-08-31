function [dpqhat_f] = dpqhatPoint2D(x0,y0,Q)
%
% The normalized differentiated Greens kernel function for free-field propagation
% Normalized: params.rho*params.omega*Q/4 = 1
%
    dpqhat_x_f = @(x,y,k) -Q/4*besselh(1, 2, k*sqrt((x - x0).^2 + (y - y0).^2))*k.*(2*x - 2*x0)./(2*sqrt((x - x0).^2 + (y - y0).^2));
    dpqhat_y_f = @(x,y,k) -Q/4*besselh(1, 2, k*sqrt((x - x0).^2 + (y - y0).^2))*k.*(2*y - 2*y0)./(2*sqrt((x - x0).^2 + (y - y0).^2));

    dpqhat_f = @(x,y,k) [dpqhat_x_f(x,y,k), dpqhat_y_f(x,y,k)];
end