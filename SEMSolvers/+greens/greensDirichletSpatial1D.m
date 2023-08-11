%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Calculate Greens Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = greensDirichletSpatial1D(Nx,xs,x0,k,tau,c)
%
% Nx,Ny: Number of wave functions in x and y dir
% M: number of grid points
% lx,ly: size x,y
% Ps: Source position
%
    if nargin == 4
        absorb_term = 0;
    else
        absorb_term = 1i*k/tau/c;
    end
    
    k_2 = k^2;
    
    if isvector(xs) && ~isrow(xs)
        xs = xs';
    end
    
    lx = xs(end);    
    P = zeros(1,length(xs));    
    V = 1;
    
    num_ignored_basis_function = 0;
    
    for nx = 0:Nx-1
        if nx == 0
            eps_x = 1;
        else
            eps_x = 2;
        end

        k_m = nx*pi/lx;
        Lambda_m = sqrt(eps_x);
        Psi_r = Lambda_m*sin(nx*pi*xs/lx);
        Psi_r0 = Lambda_m*sin(nx*pi*x0/lx)*ones(1,length(xs));

        if k-k_m == 0
            % unlimited when frequency coincides with one of
            % the natural frequencies (p. 135 Møller, Jakobsen)
            %disp('frequency coincided with natural frequency - continuing')
            num_ignored_basis_function = num_ignored_basis_function + 1;
            continue
        end
        % Greens Function
        P = P - 1/V*(Psi_r.*Psi_r0)./(k_2-k_m^2 - absorb_term);
    end
    
    P = P';
    
    fprintf('Number ignored basis function: %i\n', num_ignored_basis_function)    
end