function [i,x0] = getClosestIndex(x1d, x0_pos)
    i = find(abs(x1d - x0_pos) < eps('single'));
    if isempty(i)
        i = find(x1d > x0_pos);
        i = i(1);
    end
    x0 = x1d(i);
end