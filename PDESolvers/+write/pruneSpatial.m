function [mesh_out,p_out,dx_out] = pruneSpatial(mesh,p_in,dx_in,ppw_in,ppw_x_out)
    xprune = max(round(ppw_in/ppw_x_out),1);
    mesh_out = mesh(1:xprune:end, :);
    p_out = p_in(:,1:xprune:end,:);
    dx_out = dx_in/xprune;
end