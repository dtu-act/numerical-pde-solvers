function P_out = toStructuredGrid(P, X, Y, conn, dx, Nx, Ny, bounding_box)   
    
    x_offset = -bounding_box(1,1);
    y_offset = -bounding_box(1,2);
    
    etov = conn(:,1:3); % only extract grid points
    
    x1d = X(:);
    y1d = Y(:);
    
    P_out = zeros(Ny,Nx);
    
    for e=1:size(etov,1)
        for i=etov(e,:)
            xi = x_offset + x1d(i);
            yi = y_offset + y1d(i);
            
            i_grid = xi/dx+1;
            j_grid = yi/dx+1;
            
            % make sure the points are from a structured grid with dx
            % spacing
            %assert(abs(round(i_grid)-i_grid) < 1e-10 && abs(round(j_grid)-j_grid) < 1e-10)
            if ~(abs(round(i_grid)-i_grid) < 1e-10 && abs(round(j_grid)-j_grid) < 1e-10)
                fprintf('Point skipped: not on structured grid\n')
                continue
            end
                
            P_out(round(i_grid),round(j_grid)) = P(i);            
        end        
    end
end