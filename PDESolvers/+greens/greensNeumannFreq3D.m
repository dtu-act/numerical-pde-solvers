%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Calculate Greens Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = greensNeumannFreq3D(Nx,Ny,Nz,lx,ly,lz,Ps_all,Pr_all,tau,k,c)
    x0 = Ps_all(1,1);
    y0 = Ps_all(1,2);
    z0 = Ps_all(1,3);

    x = Pr_all(1,1);
    y = Pr_all(1,2);
    z = Pr_all(1,3);    
    
    k_2 = k.^2;
    f_m = zeros(Nx,Ny,Nz);
    
    V = lx*ly*lz;
    
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

            for nz = 0:Nz-1
                if nz == 0
                    eps_z = 1;
                else
                    eps_z = 2;
                end
                
                k_m = sqrt((nx*pi/lx).^2 + (ny*pi/ly).^2 + (nz*pi/lz).^2);
                f_m(nx+1,ny+1,nz+1) = c/2*sqrt((nx/lx).^2 + (ny/ly).^2 + (nz/lz).^2);
                Lambda_m = sqrt(eps_x*eps_y*eps_z);
                Psi_r = Lambda_m*cos(nx*pi*x/lx)*cos(ny*pi*y/ly)*cos(nz*pi*z/lz);
                Psi_r0 = Lambda_m*cos(nx*pi*x0/lx)*cos(ny*pi*y0/ly)*cos(nz*pi*z0/lz);            
                
                % Greens Function
                P = P - 1/V*(Psi_r*Psi_r0)./(k_2-k_m^2-1i*k/tau/c);
            end
        end
    end
end