function [x,y,r,s] = verticesToNodes2D(P,etov,vx,vy)
%
%    Gordon-Hall transfinite interpolation procedure [28] to map the local coordinates 
%    to the physical coordinates (x, y) using the transformation
%
%    input: P::Int,etov::Matrix{Int}, vx::Vector, vy::Vector
%
%    out: Tuple{Matrix, Matrix, Vector, Vector}
%         x,y: global x,y coordinates as matrices of size (Mp x N Matrix)
%         r,s: local coordinates on the unit triangle [-1,1]

%    See p. 80-81 "The Spectral/hp-Finite Element Method for PDE" by A. Engsig-Karup
%

    [x, y] = spectral.Nodes2D(P);   % coordinates on the equilateral triangle for polynomial of order N
    [r, s] = spectral.xytors(x, y); % map to unit triangle r,s in [-1,1] 

    % map to global coordinates for each element (triangle)
    v1 = etov(:,1);
    v2 = etov(:,2);
    v3 = etov(:,3);
    x = 0.5 * (-(r + s)*vx(v1)' + (1 + r)*vx(v2)' + (1 + s)*vx(v3)');
    y = 0.5 * (-(r + s)*vy(v1)' + (1 + r)*vy(v2)' + (1 + s)*vy(v3)');
end