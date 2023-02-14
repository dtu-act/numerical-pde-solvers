function [arr_zip] = zip2(arr1, arr2)
    assert(numel(arr1) == numel(arr1))

    arr_zip = zeros(numel(arr1()), 2);
    for i = 1:numel(arr1)
       a1 = arr1(i);
       a2 = arr2(i);
       arr_zip(i,:) = [a1,a2];
    end
end