function [P, X_orig, Y_orig] = postprocess2D(Pw, domain, acustic_params, wavenumber, scf, grid_res)

    if nargin == 5
        grid_res = 0.01;
    end
    
    lx = domain.lx;
    ly = domain.ly;
    
    % dimensions: 1xM 
    Kxws = wavenumber.Kxws; % do not use ' -> conjugate transpose
    Kyws = wavenumber.Kyws;
    Kxwr = wavenumber.Kxwr;
    Kywr = wavenumber.Kywr;

    scf_s = scf.scf_s;
    scf_r = scf.scf_r;
    
    [X,Y] = meshgrid(0:grid_res:lx, 0:grid_res:ly);
    X_orig = domain.origin(1) + X;
    Y_orig = domain.origin(2) + Y;    
    
    n_g = length(X(:));

    P = zeros(size(X));
    ones_r = ones(n_g,length(Kxwr));
    ones_s = ones(n_g,length(Kxws));
    
    phi_r = exp(-1i*( X(:).*(Kxwr.*ones_r) - (lx*(scf_r.*ones_r)).*(Kxwr.*ones_r) )).*cos( Y(:).*(Kywr.*ones_r) );
    phi_s = cos( X(:).*(Kxws.*ones_s) ).*exp( -1i*( Y(:).*(Kyws.*ones_s) - ly*(scf_s.*ones_s).*(Kyws.*ones_s) ));
    phi = [phi_r phi_s];
    
    P(1:n_g) = phi*Pw;
        
    pqhat = zeros(n_g,1);
    if domain.source.type == models.SourceType.PointSource
        pqhat(:) = domain.source.F(X_orig(:),Y_orig(:),acustic_params.k);
    end

    P(:) = P(:) + pqhat;
end