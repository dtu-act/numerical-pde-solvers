function [samples] = normalizeGRFs(samples)
    for i = 1:size(samples,1)
        if min(samples(i,:)) < -1
            samples(i,:) = samples(i,:)/abs(min(samples(i,:)));
        end
        if max(samples(i,:)) > 1
            samples(i,:) = samples(i,:)/abs(max(samples(i,:)));
        end
    end
end