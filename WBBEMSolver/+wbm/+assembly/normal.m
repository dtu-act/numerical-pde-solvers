function [nv, J] = normal(coords_pair)
    dv = coords_pair(2,:) - coords_pair(1,:);
    dv_norm = norm(dv);
    nv = [dv(2), -dv(1)]/dv_norm;
    J = dv_norm/2;   % Jacobian of boundary : dv/dksi in the case of 2D
    %fprintf("Normal vector: (%0.1f,%0.1f)\n", nv(1), nv(2))
end