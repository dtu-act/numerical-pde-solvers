function [P,XY] = greensNeumannUnstructuredSpatial2D(Nx,Ny,XY,xy0,k,tau,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Calculate Greens Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Nx,Ny: Number of wave functions in x and y dir
% M: number of grid points
% lx,ly: size x,y
% Ps: Source position
%
    if nargin == 7
        absorb_term = 1i*k/(tau*c);
    else
        absorb_term = 0;
    end
    
    X = XY(:,1);
    Y = XY(:,2);

    lx = max(X, [], 'all') - min(X, [], 'all');
    ly = max(Y, [], 'all') - min(Y, [], 'all');

    x0 = xy0(1);
    y0 = xy0(2);
    
    k_2 = k^2;  

    P = zeros(length(X), 1);
    
    V = lx*ly;
    
    num_ignored_basis_function = 0;
    
    for nx = 0:Nx-1
        if nx == 0
            eps_x = 1;
        else
            eps_x = 2;
        end

        for ny = 0:Ny-1
            if ny == 0
                eps_y = 1;
            else
                eps_y = 2;
            end

            k_m = sqrt((nx*pi/lx).^2 + (ny*pi/ly).^2);
            Lambda_m = sqrt(eps_x*eps_y);
            Psi_r = Lambda_m*cos(nx*pi*X/lx).*cos(ny*pi*Y/ly);
            Psi_r0 = Lambda_m*cos(nx*pi*x0/lx)*cos(ny*pi*y0/ly)*ones(length(X),1);

            if k-k_m == 0
                % unlimited when frequency coincides with one of
                % the natural frequencies (p. 135 Møller, Jakobsen)
                disp('frequency coincided with natural frequency - continuing')
                num_ignored_basis_function = num_ignored_basis_function + 1;
                continue
            end
            % Greens Function
            P = P - 1/V*(Psi_r.*Psi_r0)./(k_2-k_m^2 - absorb_term); %1i*1e-4
        end
    end
    
    fprintf('Number ignored basis function: %i\n', num_ignored_basis_function)
    
end