function out = BTn(dphi_x, dphi_y, coords_pair, w_vec, nw)
    [nv, J] = wbm.assembly.normal(coords_pair);

    if nargin==3
        out = J*transpose(nv(1)*dphi_x + nv(2)*dphi_y);
    else
        W = repmat(w_vec, [1,nw])';
        out = J*W.*transpose(nv(1)*dphi_x + nv(2)*dphi_y);
    end
end