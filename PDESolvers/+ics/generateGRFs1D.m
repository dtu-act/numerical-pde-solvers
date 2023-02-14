function [samples_masked, x, samples, mask, gaussian_compare] = ...
        generateGRFs1D(xminmax, sigma_0, l_0, num_of_nodes, num_of_samples, sigma0_window)
    
    if length(sigma_0) ~= length(l_0)
        raise ValueError('Size of sigma_0 and l_0 mismatch')
    end    
    
    xmin = xminmax(1);
    xmax = xminmax(2);

    x = linspace(xmin, xmax, num_of_nodes);
    dx = x(2)-x(1);
    
    [x1, x2] = meshgrid(x, x);
    distances_squared  = (x1 - x2).^2;
    covariance_matrix = zeros(num_of_nodes, num_of_nodes);
    for i = 1:length(l_0)
        % normalization factor changed from original code -> (sigma_0[i] ** 2)
        covariance_matrix = covariance_matrix + (1/(sigma_0(i)*sqrt(2*pi)) * ...
                               exp(-0.5 / (l_0(i).^2) * distances_squared));
    end
    
    mu = zeros(length(x),1);
    samples = mvnrnd(mu, covariance_matrix, num_of_samples);
    samples = reshape(samples, [], num_of_nodes);
    
%     for i = 1:size(samples,1)
%         if min(samples(i,:)) < -1
%             samples(i,:) = samples(i,:)/abs(min(samples(i,:)));
%         end
%         if max(samples(i,:)) > 1
%             samples(i,:) = samples(i,:)/abs(max(samples(i,:)));
%         end
%     end

    offset = 3.0*sigma0_window;
    x_half = 0:dx:(offset+1e-10);
    gauss_half = exp(-(x_half/sigma0_window).^2);
    mask = ones(1,size(samples,2));
    mask(1:length(x_half)) = fliplr(gauss_half);
    mask(end-length(x_half)+1:end) = gauss_half;
    
    samples_masked = samples.*mask;
    
    for i = 1:size(samples_masked,1)
        if min(samples_masked(i,:)) < -1
            samples_masked(i,:) = samples_masked(i,:)/abs(min(samples_masked(i,:)));
        end
        if max(samples_masked(i,:)) > 1
            samples_masked(i,:) = samples_masked(i,:)/abs(max(samples_masked(i,:)));
        end
    end

    gaussian_compare = exp(-(x/sigma0_window).^2);