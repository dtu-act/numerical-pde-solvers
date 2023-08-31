function out = nTB(dphi_x, dphi_y, coords_pair)
    [nv, J] = wbm.assembly.normal(coords_pair);
    out = J*(nv(1)*dphi_x + nv(2)*dphi_y);
end