function [conn,X2D,Y2D,VX,VY,etov,etoe,etof,x,y,r,s,Nk,gidx,dx] = setupMesh2D(c,fmax,P,ppw,xminmax,yminmax,do_plot)
    % Generate local mesh for reference element
    
    wavelength_min = c/fmax;
    dx = wavelength_min/ppw;
    h0 = dx*P - 1e-2; % Mesh element size (distmesh parameter)

    [VX, VY, K, etov] = setup.MeshGenDistMesh2D(xminmax,yminmax,h0,do_plot);  
    
    % mesh interconnectivity-table
    %etov = conelmtab(NeX,NeY);
    [etoe,etof] = spectral.tiConnect2D(etov,3);
    Nk = size(etov,1);
    
    [x,y] = spectral.Nodes2D(P);
    [r,s] = spectral.xytors(x,y);
    
    
    v1 = etov(:,1)';
    v2 = etov(:,2)';
    v3 = etov(:,3)';
    x = 0.5*(-(r+s)*VX(v1)+(1+r)*VX(v2)+(1+s)*VX(v3));
    y = 0.5*(-(r+s)*VY(v1)+(1+r)*VY(v2)+(1+s)*VY(v3));
    
    
    % Compute index maps for node positions
    tol = min(0.2/P^2,1e-5);
    Mp = length(s);
    % Faces
    fid1 = find( abs(s+1) < tol)';
    fid2 = find( abs(r+s) < tol)';
    fid3 = find( abs(r+1) < tol)';
    % Interior
    fint = setdiff(1:Mp,[fid1 fid2 fid3]);
    % Reorder local nodes to match the chosen ordering in terms of vertex nodes,
    % edge nodes and interior nodes.
    Mpf = P+1;
    LocalReorder = [1 Mpf Mp fid1(2:Mpf-1) fid2(2:Mpf-1) fid3(Mpf-1:-1:2) fint];
    if length(LocalReorder)~=length(r(:))
       error('tol is probably not correct');
    end
        
    x = x(LocalReorder,:);
    y = y(LocalReorder,:);
    r = r(LocalReorder);
    s = s(LocalReorder);
    
    %% GENERATE CONNECTIVITY TABLE
    % create connectivity list with global numbers for each local degree of
    % freedom
    [conn,gidx,X2D,Y2D] = meshing.ConstructConnectivityTable2D(VX,VY,etov,P,x,y,etoe,etof);
    %[conn, X2D, Y2D, gidx] = meshing.connectivityTable2D(P,x,y,VX,VY,etov,etoe,etof);
end
