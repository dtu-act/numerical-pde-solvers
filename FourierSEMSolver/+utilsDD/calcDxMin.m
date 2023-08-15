function [dx] = calcDxMin(x1d)
    % assume x1d is sorted
    dx = min(x1d(2:end) - x1d(1:end-1));
end