function [A, F] = assemblyInterface2D(domain, wavenumber, scf, wavecount, trunc_par, params)

    import models.*
    import wbm.assembly.*
    
    rho = params.rho;
    omega = params.omega;
    
    index_to_coords2D = @(indexes, coords) coords(indexes(:)',:);    
    coords_bound = index_to_coords2D(domain.nodes_interface, domain.gcoord);
    
    nw = wavecount.nw;
    A = zeros(nw, nw);
    F = zeros(nw,1);

    nbounds = size(domain.nodes_interface,1);
    
    % MWR on boundary using Gaussian quadrature
    for i=0:nbounds-1
        coords_pair = coords_bound(2*i+1:2*i+2,:); % -> [1:2,:], [3:4,:], [5:6,:], ...        
        %disp(coords_pair)
        
        [phi, dphi_x, dphi_y, XY, ~, w_vec] = ...
            wbm.wavebasis2D(domain, coords_pair, wavenumber, scf, wavecount, trunc_par, params.k);

        % velocity
        dpqhat = wbm.sourceAtBoundary(i,XY,params.k,domain,models.BoundaryCondition.Velocity);        
        A(:) = A + 1i/(rho*omega)*PhiT(phi,w_vec,nw) * nTB(dphi_x, dphi_y, coords_pair);
        F(:) = F - 1i/(rho*omega)*PhiT(phi,w_vec,nw) * nTdPq(dpqhat(:,1),dpqhat(:,2),coords_pair);
        
        % pressure
        pqhat = wbm.sourceAtBoundary(i,XY,params.k,domain,models.BoundaryCondition.Pressure);        
        A(:) = A - 1i/(rho*omega)*BTn(dphi_x,dphi_y,coords_pair) * Phi(phi,w_vec,nw);
        F(:) = F + 1i/(rho*omega)*BTn(dphi_x,dphi_y,coords_pair,w_vec,nw) * pqhat;
    end
end