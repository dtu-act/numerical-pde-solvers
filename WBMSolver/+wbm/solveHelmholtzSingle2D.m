function [P, XY, wavecount] = solveHelmholtzSingle2D(domain, params, wx_max, wy_max, grid_res)
    
    if nargin == 4
        grid_res = 0.1;
    end
    
    [wavenumber, scf, wavecount] = wbm.wavenumber_gen2D(domain, params, wx_max, wy_max);

    trunc_par_x = ceil(wx_max * pi / (domain.lx*params.k));
    %fprintf("trunc_par_x: %i\n", trunc_par_x)

    [A, F] = wbm.assembly.assembly2D(domain, wavenumber, scf, wavecount, trunc_par_x, params);

    Pw = A \ F;
    
    [P, xs, ys] = wbm.postprocess2D(Pw, domain, params, wavenumber, scf, grid_res);
    XY = [xs(:),ys(:)];
end

% function run2D()::Tuple{Matrix, Vector, Vector}
%     params = AcousticModels.AcousticParameters(500,340,1.225)
%     trunc_par = 3      # Wave function truncation parameter
% 
%     l_x = 0.5
%     l_y = 0.5
%     gcoord=[0.5 0.125; 1.0 0.0; 1 0.5; 0.5 0.375] # coordinates at boundaries
%     origin = [0.5, 0]
% 
%     bound_nodes = [1 2; 2 3; 3 4; 4 1]
%     bound_velocities = [ 0 0; 1 0; 0 0; 0 0]
% 
%     domain = AcousticModels.DomainWBM2D(l_x, l_y, gcoord, origin, bound_nodes, AcousticModels.velocity_t, bound_velocities)
%     
%     wavenumber, scf, wavecount = WBUtils.wavenumber_gen2D(domain, trunc_par, params)
% 
%     (A, F) = MatrixAssembly.assemblyWBM2D(domain, wavenumber, scf, wavecount, trunc_par, params)
%     
%     print("A size: $(size(A))\n")
% 
%     Pw = A \ F
% 
%     P, xs, ys = MatrixAssembly.postprocessWBM2D(domain, params, wavenumber, scf, wavecount, Pw, grid_res=0.001)
% 
%     return P, origin[1].+xs, origin[2].+ys
% end