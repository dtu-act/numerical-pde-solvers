function [C, Fij] = assemblyCoupling2D(i,j,domains,wavenumbers,scfs,wavecounts,trunc_par,params_acc)
    import models.*
    import wbm.assembly.*

    rho = params_acc.rho;
    omega = params_acc.omega;

    domain_i = domains{i}; domain_j = domains{j};
    wavenumber_i = wavenumbers{i}; scf_i = scfs{i}; wavecount_i = wavecounts{i};
    wavenumber_j = wavenumbers{j}; scf_j = scfs{j}; wavecount_j = wavecounts{j};
    nw_i = wavecount_i.nw; nw_j = wavecount_j.nw;
    
    index_to_coords2D = @(indexes, coords) coords(indexes(:),:);    
    coords_interface_i = index_to_coords2D(domain_i.nodes_interface', domain_i.gcoord);
    coords_interface_j = index_to_coords2D(domain_j.nodes_interface', domain_j.gcoord);

    C = zeros(nw_i, nw_j);
    Fij = zeros(nw_i,1);

    assert(length(domain_i.nodes_interface) == length(domain_j.nodes_interface))
    
    nbounds = size(domain_i.nodes_interface,1);
    
    % MWR on boundary using Gaussian quadrature
    for k=0:nbounds-1
        coords_pair_i = coords_interface_i(2*k+1:2*k+2,:); % -> [1:2,:], [3:4,:], [5:6,:], ...
        coords_pair_j = coords_interface_j(2*k+1:2*k+2,:); % -> [1:2,:], [3:4,:], [5:6,:], ...
        
%         disp('alpha')
%         disp(coords_pair_i)
%         disp('beta')
%         disp(coords_pair_j)
        
        [phi_i, dphi_x_i, dphi_y_i, ~, ~, w_vec_i] = ...
            wbm.wavebasis2D(domain_i, coords_pair_i, wavenumber_i, scf_i, wavecount_i, trunc_par, params_acc.k);

        [phi_j, dphi_x_j, dphi_y_j, XY_j, ~, w_vec_j] = ...
            wbm.wavebasis2D(domain_j, coords_pair_j, wavenumber_j, scf_j, wavecount_j, trunc_par, params_acc.k);
              
        % velocity
        dpqhat_j = wbm.sourceAtBoundary(k,XY_j,params_acc.k,domain_j,models.BoundaryCondition.Velocity);
        C(:) = C + 1i/(rho*omega)*PhiT(phi_i,w_vec_i,nw_i) * nTB(dphi_x_j, dphi_y_j, coords_pair_j);
        Fij(:) = Fij - 1i/(rho*omega)*PhiT(phi_i,w_vec_i,nw_i) * nTdPq(dpqhat_j(:,1),dpqhat_j(:,2),coords_pair_j);
        
        % pressure
        pqhat_j = wbm.sourceAtBoundary(k,XY_j,params_acc.k,domain_j,models.BoundaryCondition.Pressure);
        C(:) = C + 1i/(rho*omega)*BTn(dphi_x_i,dphi_y_i,coords_pair_i) * Phi(phi_j,w_vec_j,nw_j);
        Fij(:) = Fij - 1i/(rho*omega)*BTn(dphi_x_i,dphi_y_i,coords_pair_i,w_vec_i,nw_i) * pqhat_j;
    end
end