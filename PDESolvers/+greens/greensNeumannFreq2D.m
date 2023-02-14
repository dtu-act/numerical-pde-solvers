%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Calculate Greens Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, f_m] = greensNeumannFreq2D(Nx,Ny,lx,ly,Ps,Pr,k,c,tau)
    
    if nargin == 8
        absorb_term = 0;
    else
        absorb_term = 1i*k/tau/c;
    end
    
    x0 = Ps(1); y0 = Ps(2);
    x = Pr(1); y = Pr(2);   
    
    k_2 = k.^2;
    f_m = zeros(Nx,Ny);
    
    V = lx*ly;
    
    P = 0;
    
    % Loop over all combinations of modes
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
            f_m(nx+1,ny+1) = c/2*sqrt((nx/lx).^2 + (ny/ly).^2);
            Lambda_m = sqrt(eps_x*eps_y);
            Psi_r = Lambda_m*cos(nx*pi*x/lx)*cos(ny*pi*y/ly);
            Psi_r0 = Lambda_m*cos(nx*pi*x0/lx)*cos(ny*pi*y0/ly);

            % Greens Function
            P = P - 1/V*(Psi_r*Psi_r0)./(k_2-k_m^2 - absorb_term);
        end
    end
end