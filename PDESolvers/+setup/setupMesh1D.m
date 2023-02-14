function [vx,etov,conn,rs,x1d] = setupMesh1D(c,fmax,Porder,ppw,xminmax)
    xl = xminmax(1);
    xu = xminmax(2);
    
    % create mesh (no higher order nodes)
    [vx, etov] = sem.setupUniformGrid1D(c,fmax,Porder,ppw,xl,xu);
    
    % local nodes as Legendre-Gauss-Lobatto (LGL) (p. 63)
    rs = spectral.JacobiGL(0, 0, Porder);
    
    %(nodes, etov) = reorderUniformGrid1D(vx, etov) % not needed in this case, but in general mesh1D requires sorted input
    conn = meshing.constructConnectivity1D(etov, Porder);
    
    % insert additional elements to the mesh, so that the LGL nodes gets inserted
    % into each local element
    x1d =  meshing.mesh1D(vx, rs);