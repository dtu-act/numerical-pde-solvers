function L = Lagrange1D(rs, xs)
%
% Input: 
%   rs: absissas, e.g. Gauss-Legendre nodal points
%   xs: the grid point to interpolate
% Output: 
%   L: matrix of the interpolating Lagrange function L_i(x_j) = psi_i(x_j)
%   with rows corresponding to polynomial orders 0,1,2,..N and
%   columns corresponding to the grid index
%
    N = length(rs)-1;
    V = spectral.Vandermonde1D(N, rs);
    PSI = spectral.Vandermonde1D(N, xs);
        
    L = zeros(N+1,length(xs));
    L(:) = V'\PSI';
    
    %PSI_absisas = spectral.Vandermonde1D(N, rs);
    %L_abs = V'\PSI_absisas';
    %plot(xs, L(5,:),'-'); hold on; plot(rs, L_abs(5,:),'o'); hold off
    
%     for n=0:N
%         Vinv = inv(V');
%         L(:) = L(:) + Vinv(:,n+1).*spectral.JacobiP(xs,0,0,n);
%     end
end