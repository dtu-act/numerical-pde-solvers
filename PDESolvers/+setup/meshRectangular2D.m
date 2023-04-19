function [VX, VY, K, EToV] = meshRectangular2D(xminmax,yminmax,h0,do_plot)
    % function [VX, VY, K, EToV] = MeshGenDistMesh2D()
    % Purpose : Generate 2D square mesh using DistMesh;
    % By Allan P. Engsig-Karup
    % Modified by Finnur Pind
    
    % Parameters to set/define
    % fd Distance function for mesh boundary
    % fh Weighting function for distributing elements
    % h0 Characteristic length of elements
    % Bbox Bounding box for mesh
    % param Parameters to be used in function call with DistMesh
    
    xl = xminmax(1);
    xu = xminmax(2);
    
    yl = yminmax(1);
    yu = yminmax(2);
    
    cornerp = [xl,yl;xu,yl;xl,yu;xu,yu];
    fd = @(p) drectangle(p,xl,xu,yl,yu);
    fh = @huniform;
    Bbox = [xl yl; xu yu];
    
    % Call distmesh
    [Vert,EToV] = distmesh2d(fd,fh,h0,Bbox,cornerp,do_plot);
    VX = Vert(:,1)'; VY = Vert(:,2)';
    Nv = length(VX); K = size(EToV,1);
    
    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));
    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);
end