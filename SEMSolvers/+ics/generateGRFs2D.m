function [samples_masked] = generateGRFs2D(sigma_0, l_0, XX, YY, num_samples, sigma0_win)
    if length(sigma_0) ~= length(l_0)
        raise ValueError('Size of sigma_0 and l_0 mismatch')
    end
    
    num_nodes_x = size(XX,1);
    num_nodes_y = size(YY,2);
    num_of_nodes = num_nodes_x * num_nodes_y;
    dx = XX(1,2) - XX(1,1);
    dy = YY(2,1) - YY(1,1);

    x = XX(:);
    y = YY(:);
    [x1, x2] = meshgrid(x, x); %, indexing='ij'
    [y1, y2] = meshgrid(y, y); % , indexing='ij'
    distances_squared  = ((x1 - x2).^2 + (y1 - y2).^2);
    covariance_matrix = zeros(num_of_nodes, num_of_nodes);
    for i = 1:length(l_0)
        % normalization factor changed from original code: 
        % From (sigma_0[i] ** 2) to (1/(sigma_0[mode_index]*np.sqrt(2*np.pi))
        covariance_matrix = covariance_matrix + (1/(sigma_0(i)*sqrt(2*pi)) * ...
                               exp(- 0.5 / (l_0(i).^2) * distances_squared));
    end

    mu = zeros(length(x),1);
    fprintf('Calculating %i GRF samples...\n', num_samples)
    samples = mvnrnd(mu, covariance_matrix, num_samples);
    samples = reshape(samples, [], num_nodes_x, num_nodes_y);
    %samples = samples.reshape([-1, num_of_nodes_x, num_of_nodes_y])

    % normalize
    for i = 1:size(samples,1)
        if min(samples(i,:)) < -1
            samples(i,:) = samples(i,:)/abs(min(samples(i,:)));
        end
        if max(samples(i,:)) > 1
            samples(i,:) = samples(i,:)/abs(max(samples(i,:)));
        end
    end

    offset = sigma0_win*3;
    x_half = 0:dx:offset+1e-10;
    y_half = 0:dy:offset+1e-10;
    gauss_half_x = exp(-(x_half/sigma0_win).^2);
    gauss_half_y = exp(-(y_half/sigma0_win).^2);
    
    gauss_xbound = repmat(gauss_half_x, num_nodes_y, 1);
    gauss_ybound = repmat(gauss_half_y', 1, num_nodes_x);

    mask = ones(size(samples,2), size(samples,3));
    
    Nx = length(x_half);
    Ny = length(y_half);
    mask(:,1:Nx) = mask(:,1:Nx).*flip(gauss_xbound,2);
    mask(:,end-Nx+1:end) = mask(:,end-Nx+1:end).*gauss_xbound;
    mask(1:Ny,:) = mask(1:Ny,:).*flip(gauss_ybound,1);
    mask(end-Ny+1:end,:) = mask(end-Ny+1:end,:).*gauss_ybound;
    
    mask_ = ones(size(samples));
    for i=1:num_samples
        mask_(i,:,:) = mask;
    end
    samples_masked = samples.*mask_;
    %samples_masked = reshape(samples_masked, size(samples_masked,1), []);