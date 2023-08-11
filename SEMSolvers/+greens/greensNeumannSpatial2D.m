function [P, X, Y] = greensNeumannSpatial2D(Nx,Ny,lx,ly,Ne,xy0,k,tau,c)
    
    if nargin == 7
        absorb_term = 0;
    else
        absorb_term = 1i*k/tau/c;
    end
    
    x0 = xy0(1);
    y0 = xy0(2);
    
    k_2 = k^2;
    
    xs = linspace(0,lx,Ne+1);
    ys = linspace(0,ly,Ne+1);
    
    [X, Y] = meshgrid(xs,ys);
    P = zeros(length(ys), length(xs));
    
    V = lx*ly;
    
    num_ignored_basis_function = 0;
    % Loop over all combinations of modes
    for i=1:length(xs)
        for j=1:length(ys)     
            x = xs(i);
            y = ys(j);
            
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
                    Psi_r  = Lambda_m*cos(nx*pi*x/lx)*cos(ny*pi*y/ly);
                    Psi_r0 = Lambda_m*cos(nx*pi*x0/lx)*cos(ny*pi*y0/ly);

                    if k-k_m == 0
                        % unlimited when frequency coincides with one of
                        % the natural frequencies (p. 135 Møller, Jakobsen)
                        num_ignored_basis_function = num_ignored_basis_function + 1;
                        continue
                    end
                    % Greens Function
                    P(i,j) = P(i,j) - 1/V*(Psi_r*Psi_r0)./(k_2-k_m^2 - absorb_term);
                end
            end
        end
    end
    
    fprintf('Number ignored basis function: %i\n', num_ignored_basis_function)
    
end