function out = nTdPq(pqhat_x,pqhat_y,coords_pair)
    [nv, J] = wbm.assembly.normal(coords_pair);
    out = J*(nv(1)*pqhat_x + nv(2)*pqhat_y);    
end