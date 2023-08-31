function out = PhiT(phi,w_vec,nw)
    W = repmat(w_vec, [1,nw]);
    out = transpose(W.*phi);
end