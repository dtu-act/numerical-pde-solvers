function [VX, VY, EToV, dx] = meshPolygon2D(c,fmax,P,ppw,pv,do_plot)
    % Purpose : Generate 2D L-shape mesh using DistMesh;
    
    % Parameters to set/define
    % fd Distance function for mesh boundary
    % fh Weighting function for distributing elements
    % h0 Characteristic length of elements
    % Bbox Bounding box for mesh
    % param Parameters to be used in function call with DistMesh
    
    wavelength_min = c/fmax;
    dx = wavelength_min/ppw;
    h0 = dx*P - 1e-2; % Mesh element size (distmesh parameter)
    
    xminmax = [min(pv(:,1)),max(pv(:,1))];
    yminmax = [min(pv(:,2)),max(pv(:,2))];

    xl = xminmax(1);
    xu = xminmax(2);
    yl = yminmax(1);
    yu = yminmax(2);

    bbox = [xl yl; xu yu];
    if size(pv,1) == 4
        fd = @(p) drectangle(p,xl,xu,yl,yu);
    else
        fd = { 'l_dpolygon', [], pv };
    end    
    fh = @(p) ones(size(p,1),1);

    % Call distmesh
    [Vert,EToV] = distmesh( fd, fh, h0, bbox, pv );
    VX = Vert(:,1)'; VY = Vert(:,2)';
    if do_plot
        patch( 'vertices', Vert, 'faces', EToV, 'facecolor', [.9, .9, .9] )
    end

    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));
    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);
end