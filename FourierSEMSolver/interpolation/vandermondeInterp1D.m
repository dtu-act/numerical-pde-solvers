function [poly_interp_f] = vandermondeInterp1D(x,y)
    num_interp = length(x);
    
    V = [];
    for j=1:num_interp
        V = [V; x(j).^(0:num_interp-1)];
    end
    %V = spectral.Vandermonde1D(num_interp-1,x);
    
    weights = V \ y';
    
    function [z] = poly_interp(x)
        z = zeros(length(x),1);
        for i = 1:length(x)
            z(i) = x(i).^(0:num_interp-1)*weights;
        end
    end

    poly_interp_f = @poly_interp;
end