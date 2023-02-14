function [etov] = setupEToV1D(vx)
    etov = zeros(length(vx)-1, 2);
    for i=1:length(vx)-1
        etov(i,:) = [i, i+1];
    end
end