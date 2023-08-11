function [phi, dphi_x, dphi_y, Xi_origin,r,W] = wavebasis2D(domain, bcoord_pair, wavenumber, scf, wavecount, trunc_par, k)
%
%    Construct the wave functions and its derivatives evaluated at the Gaussian quadrature points
%

    lx = domain.lx;
    ly = domain.ly;
    origin = domain.origin;
    
    % dimensions: 1xM    
    Kxws = wavenumber.Kxws;
    Kyws = wavenumber.Kyws;
    Kxwr = wavenumber.Kxwr;
    Kywr = wavenumber.Kywr;
    
    if size(Kxws,1) ~= 1 || size(Kyws,1) ~= 1 || size(Kxwr,1) ~= 1 || size(Kywr,1) ~= 1
        error("Wrong wave number dimension")
    end
    
    scf_s = scf.scf_s;
    scf_r = scf.scf_r;

    % dimensions: scalar
    nwr = wavecount.nwr;
    nws = wavecount.nws;

    % dimensions: 1x2
    dv = bcoord_pair(2,:) - bcoord_pair(1,:);

    % dimensions: scalar
    dv_norm = norm(dv);
    
    %  Gaussian quadrature at each boundary 
    n_quadp = max(2,ceil(6*trunc_par*dv_norm*k/pi));
    [r,W] = spectral.JacobiGQ(0,0,n_quadp-1);
    
    Xi = zeros(n_quadp,2);
    Xi_origin = zeros(n_quadp,2);
    for n=1:n_quadp
        H = [ (1-r(n)) (1+r(n)) ]./2;
        xi = H*bcoord_pair;
        xi_offset = xi - origin;
        Xi(n,1:2) = xi_offset;
        Xi_origin(n,1:2) = xi;
    end

    % ============ Construction of wave function ==============
    phi_r   = zeros(n_quadp,nwr);
    phi_s   = zeros(n_quadp,nws);
    dphi_rx = zeros(n_quadp,nwr);
    dphi_ry = zeros(n_quadp,nwr);
    dphi_sx = zeros(n_quadp,nws);
    dphi_sy = zeros(n_quadp,nws);
    
    %fprintf("n_quadp %i\n", n_quadp)
    %fprintf("nw %i\n", nwr + nws)
    
    % NOTE: only x,y pairwise (since on a 1D boundary): [(x_0,y_0), (x_1,y_1), ...]
    
    ones_col = ones(n_quadp,1);
    
    phi_s(:) = cos( Xi(:,1)*Kxws ).* ...
        exp(-1i*( Xi(:,2)*Kyws - (ly*ones_col*scf_s).*(ones_col*Kyws) ));
    phi_r(:) = exp(-1i*( Xi(:,1)*Kxwr - (lx*ones_col*scf_r).*(ones_col*Kxwr) )).* ...
        cos( Xi(:,2)*Kywr );
    
    dphi_sx(:) = (-ones_col*Kxws).* ...
        sin( Xi(:,1)*Kxws ).* ...
        exp(-1i*( Xi(:,2)*Kyws - (ly*ones_col*scf_s).*(ones_col*Kyws) ));
    dphi_sy(:) = (-1i*ones_col*Kyws).* ...
        cos( Xi(:,1)*Kxws ).* ...
        exp(-1i*( Xi(:,2)*Kyws - (ly*ones_col*scf_s).*(ones_col*Kyws) ));
    
    dphi_rx(:) = (-1i*ones_col*Kxwr).* ...
        exp(-1i*( Xi(:,1)*Kxwr - (lx*ones_col*scf_r).*(ones_col*Kxwr) )).* ...
        cos( Xi(:,2)*Kywr );
    dphi_ry(:) = (-ones_col*Kywr).* ...
        exp(-1i*( Xi(:,1)*Kxwr - (lx*ones_col*scf_r).*(ones_col*Kxwr) )).* ...
        sin( Xi(:,2)*Kywr );
    
    phi    = [phi_r phi_s];
    dphi_x = [dphi_rx dphi_sx];
    dphi_y = [dphi_ry dphi_sy];
end